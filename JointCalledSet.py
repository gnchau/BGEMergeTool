import json
import pandas as pd
import hail as hl
from datetime import date
import numpy as np
from bokeh.plotting import figure, show
from bokeh.io import output_notebook
from hail.fs.hadoop_fs import hadoop_open


class JointCalledSet:
    """
    Represents a joint set of genetic variant calls obtained from exome sequencing (exo)
    and imputation (imp) datasets. Provides methods for loading, merging, and analyzing these datasets.

    Args:
        exo_path (str): Path to the exome variant call file in VCF format.
        imp_path (str): Path to the imputation variant call file in VCF format.
        output_directory (str): Directory where the output files will be saved.
        priority (str, optional): The priority dataset ('exome' or 'imputation'). Defaults to 'exome'.

    Attributes:
        exo (hl.MatrixTable): Hail MatrixTable representing the exome dataset.
        imp (hl.MatrixTable): Hail MatrixTable representing the imputation dataset.
        out_path_prefix (str): Prefix for the output file paths.
        priority (str): The priority dataset ('exome' or 'imputation').
        merged (hl.MatrixTable): Merged dataset with priority given to the specified dataset.
        overlap (hl.MatrixTable): Dataset representing the overlap of variants and samples between exome and imputation datasets.
    """
    def __init__(self, 
                 exo_path, 
                 imp_path, 
                 output_directory, 
                 log_outname=f'merging-tool-{date.today().isoformat()}-{np.random.randint(100,999)}.log',
                 priority='exome',
                 keep_all_fields=False,
                 impute_nonsimilar_samples=False):
        self.exo_path, self.imp_path = exo_path, imp_path
        self.load_sets(exo_path, imp_path)
        self.log = ''
        _msg = ''
                
        if output_directory[-1] != '/':
            raise ValueError('Output Directory must be a directory!')
        else:
            self.out_path_prefix = output_directory

        self.log_outpath = f'{self.out_path_prefix}{log_outname}'
        print(f'Logging session at {self.log_outpath}.')

        self.exo = self.exo.annotate_entries(
            GT = hl.vds.lgt_to_gt(self.exo.LGT, self.exo.LA)
        )

        self.priority = priority
        
        if not keep_all_fields:
            # drop these but we might want them later
            fields_to_remove_exo = list(self.exo.entry.dtype.fields)
            fields_to_remove_exo.append('info')
            fields_to_remove_exo.remove('GT')

            fields_to_remove_imp = list(self.imp.entry.dtype.fields)
            fields_to_remove_imp.append('info')
            fields_to_remove_imp.remove('GT')
            self.exo = self.exo.drop(*fields_to_remove_exo)
            self.imp = self.imp.drop(*fields_to_remove_imp)
        
        if impute_nonsimilar_samples:
            samples_in_imp = set(self.imp.s.collect())
            samples_in_exo = set(self.exo.s.collect())
            non_intersecting_samples = samples_in_imp.symmetric_difference(samples_in_exo)

            for sample in non_intersecting_samples:
                self.imp = self.imp.annotate_entries(**{sample: hl.missing(hl.tcall)})

            joined_dataset = self.imp.union_rows(self.exo)

            imputed_gt_value = hl.call(0, 0)
            joined_dataset = joined_dataset.annotate_entries(**{sample: hl.coalesce(joined_dataset[sample], imputed_gt_value) for sample in non_intersecting_samples})

        self._num_exo_sites, self._num_exo_samples = self.exo.count()
        self._num_imp_sites, self._num_imp_samples = self.imp.count()
        
        _msg += f'Loaded {self._num_exo_sites} sites and {self._num_exo_samples} samples from {self.exo_path}.\n'
        _msg += f'Loaded {self._num_imp_sites} sites and {self._num_imp_samples} samples from {self.imp_path}.\n'
        print(_msg)
        
        if self._num_exo_samples != self._num_imp_samples:
            samples_diff = list(set(self.exo.s.collect()) - set(self.imp.s.collect()))
            diff_path = f'{self.out_path_prefix}no_sample_match.txt'
            self.export_flat(samples_diff, diff_path)
            _msg += f'Warning: {len(samples_diff)} have no match! Saving to {diff_path}.\n'
            print(f'Warning: {len(samples_diff)} have no match! Saving to {diff_path}.')
        else:
            _msg += 'Success! All samples matched.\n'
            print('Success! All samples matched.')
                                   
        # takes union of rows (sites), priority on exomes by default.
        _msg += 'Merging...\n'
        print('Merging...')
        
        self.merge()
        
        #self.merged = self.merged.repartition(1000)
        _msg += f'Done merging! Priority has been given to {self.exo_path}.\nTherefore {self.exo_path} variant calls will be preserved and {self.imp_path} \
        will be overwritten.\n'
        print(f'Done merging! Priority has been given to {self.exo_path}.\nTherefore {self.exo_path} variant calls will be preserved and {self.imp_path} will be overwritten.')

        # also calculates overlap matrix: overlap on both variants and samples
        self.calc_overlap(output_name=None)
        
        self.update_required = True
        self.exome_split = None
        self._multi_shape = None
        self.PRIORITIES = {'exome', 'imputation'}
        self.OUT_TYPES = {'vcf', 'PLINK', 'plink', 'hail'}
        self._shape = None
        self.samples = None
        self.sample_conc = None
        self.site_conc = None
        self.concordance_results = None

        _msg += f'Merged for a total of {self.shape[0]} sites and {self.shape[1]} samples.'
        print(f'Merged for a total of {self.shape[0]} sites and {self.shape[1]} samples.')
        
        self.log += _msg
    
    def compute_nonref_concordance(self, matrix):
        """
        Compute non-reference concordance from a concordance matrix.
        
        Args:
            matrix (list): A 5x5 concordance matrix.
            
        Returns:
            float: The non-reference concordance rate.
        """
        submatrix = [row[2:] for row in matrix[2:]]
        b, c = submatrix[0][1], submatrix[0][2]
        d, e, f = submatrix[1][0], submatrix[1][1], submatrix[1][2]
        g, h, i = submatrix[2][0], submatrix[2][1], submatrix[2][2]
        numerator = e + i
        denominator = (c + g) * 2 + b + d + e + f + h + i
        return numerator / denominator if denominator > 0 else 0

    def compute_concordance_metrics(self, matrix):
        """
        Compute various concordance metrics from a concordance matrix.
        
        Args:
            matrix (list): A 5x5 concordance matrix.
            
        Returns:
            dict: A dictionary containing concordance metrics.
        """
        genotype_rows = matrix[2:5] 
        n_concordant = sum(matrix[i][i] for i in range(2, 5))
        
        total_obs = sum(sum(row[2:5]) for row in genotype_rows)
        n_discordant = total_obs - n_concordant
        concordance_rate = n_concordant / total_obs if total_obs > 0 else 0
        
        nonref_concordance = self.compute_nonref_concordance(matrix)
        het_to_ref = matrix[3][2]  
        het_to_het = matrix[3][3] 
        
        return {
            'n_discordant': n_discordant,
            'n_total': total_obs,
            'concordance_rate': concordance_rate,
            'nonref_concordance': nonref_concordance,
            'het_to_ref': het_to_ref,
            'het_to_het': het_to_het
        }
        
    def load_sets(self, exo_path, imp_path):
        """
        Loads exome and imputation datasets into hail format.
        Args:
            exo_path (str): Path to the exome variant call file in VCF format.
            imp_path (str): Path to the imputation variant call file in VCF format.
        """
        if exo_path[-4:] == '.bgz' or exo_path[-3:] == '.gz':
            self.exo = hl.import_vcf(exo_path, 
                           force_bgz=True,
                           call_fields=['LGT'])
        else:
            self.exo = hl.import_vcf(exo_path, 
                           call_fields=['LGT'])
        
        if imp_path[-4:] == '.bgz' or imp_path[-3:] == '.gz':
            self.imp = hl.import_vcf(imp_path, 
                               force_bgz=True,)
        else:
            self.imp = hl.import_vcf(imp_path)
    
    def set_priority(self, new):
        """
        sets a new priority for the dataset
        Args:
            new (str): New priority ('exome' or 'imputation').
        """
        if not new in self.PRIORITIES:
            raise ValueError(f'New priority must be one of {self.PRIORITIES}.')
        self.priority = new
    
    def describe(self):
        """
        Calls describe on the merged mt.
        """
        self.merged.describe()
        
    def export_table(self, out_type, file_name_prefix='merged', overwrite=False):
        """
        exports the merged dataset in the specified format.
        Args:
            out_type (str): Output format ('vcf', 'PLINK', 'plink', 'hail').
            overwrite (bool, optional): If True, overwrite existing output files.
        """
        if not out_type in self.OUT_TYPES:
            raise ValueError(f'Export type must be in {self.OUT_TYPES}.')
        if out_type == 'hail':
            export_outpath = f'{self.out_path_prefix}{file_name_prefix}.mt'
            self.merged.write(export_outpath)
        elif out_type == 'vcf':
            export_outpath = f'{self.out_path_prefix}{file_name_prefix}.vcf.bgz'
            hl.export_vcf(self.merged, export_outpath)
        else:
            export_outpath = f'{self.out_path_prefix}{file_name_prefix}'
            hl.export_plink(self.merged, export_outpath, ind_id=self.merged.s)
            sites = self.generate_sites_file(export=True)
            
        _msg = f'Exported merged set to {export_outpath} with {self.shape[0]} sites and {self.shape[1]} samples.'
        print(_msg)
        
        self.log += _msg + '\n'
            
    def export_flat(self, itm, fname):
        """
        exports a list or string to a flat file. NOTE: requires GCS connector.
        Args:
            itm (obj): Item to be exported.
            fname (str): File name for the exported list.
        """
        if isinstance(itm, list):
            pd.DataFrame(itm).to_csv(fname, sep='\t', index=False, header=False)
        else:
            with hadoop_open(fname, 'w') as f: 
                f.write(itm)
            
    def export_samples(self, sample_file_prefix='jc_samples'):
        """
        exports the list of matching samples to a file.
        Args:
            sample_file_prefix (str, optional): Prefix for the sample file name.
        """
        if not self.samples: 
            self.samples = self.merged.s.collect()
        sample_out_path = f'{self.out_path_prefix}{sample_file_prefix}.txt'
        self.export_flat(self.samples, sample_out_path)
        
        _msg = f'Matching sample list saved to {sample_out_path}.'
        print(_msg)
        self.log += _msg + '\n'
            
    def split_multi(self):
        """
        splits multi-allelic variants and updates the internal dataset [self.exome_split].
        """
        if not self.exome_split:
            bi = self.exo.filter_rows(hl.len(self.exo.alleles) == 2)
            bi = bi.annotate_rows(a_index=1, was_split=False, old_locus=bi.locus, old_alleles=bi.alleles)
            multi = self.exo.filter_rows(hl.len(self.exo.alleles) > 2)
            split = hl.split_multi(multi)
            split = split.union_rows(bi).drop(*['a_index', 'was_split', 'old_locus', 'old_alleles'])
            self.exome_split = self.merged.union_rows(split)
        else:
            _msg = 'JointCalledSet has already been split!'
            print(_msg)
            self.log += _msg + '\n'
            
        if not self._multi_shape:
            self.count_multi_split()
            
        _msg = f'Before multi-allelic split, exome has {self._num_exo_sites} variants and {self._num_exo_samples} samples.'
        print(_msg)
        self.log += _msg + '\n'
        
        _msg = f'After multi-allelic split, exome has {self._multi_shape[0]} variants and {self._multi_shape[1]} samples.'
        print(_msg)
        self.log += _msg + '\n'
        
    def count_multi_split(self):
        self._multi_shape = self.exome_split.count()

    def concordance(self, output_name=None, filtered=False, output_json=True, graph_out=False):
        """
        Calculates concordance between exome and imputation datasets, with detailed analysis across 
        allele count (AC) and minor allele frequency (MAF) bins.
        
        Args:
            output_name (str or list, optional): Output file name or list of output file names.
            filtered (bool, optional): If True, assumes datasets are already filtered. Defaults to False.
            output_json (bool, optional): If True, outputs results to a JSON file. Defaults to True.
            graph_out (bool, optional): If True, generates histograms. Defaults to False.
        
        Returns:
            dict: A dictionary containing concordance results.
        """
        _msg = 'Calculating concordance...\n'
        print(_msg)
        self.log += _msg + '\n'
        
        if not filtered:
            # Filter to bi-allelic sites
            exome_filtered = self.exo.filter_rows(hl.len(self.exo.alleles) == 2)
            imputation_filtered = self.imp.filter_rows(hl.len(self.imp.alleles) == 2)
            
            n_samples = exome_filtered.count_cols()
            _msg = f"Total samples: {n_samples}"
            print(_msg)
            self.log += _msg + '\n'
            
            # Annotate with AC and MAF
            exome_filtered = exome_filtered.annotate_rows(
                AC=hl.agg.sum(exome_filtered.GT.n_alt_alleles()),
                MAF=hl.agg.sum(exome_filtered.GT.n_alt_alleles()) / (2 * n_samples)
            )
            
            imputation_filtered = imputation_filtered.annotate_rows(
                AC=hl.agg.sum(imputation_filtered.GT.n_alt_alleles()),
                MAF=hl.agg.sum(imputation_filtered.GT.n_alt_alleles()) / (2 * n_samples)
            )
            
            # Find intersection of variants
            exome_variants = exome_filtered.rows().select()._key_by_assert_sorted('locus', 'alleles')
            imputation_variants = imputation_filtered.rows().select()._key_by_assert_sorted('locus', 'alleles')
            
            intersection = exome_variants.semi_join(imputation_variants)
            
            intersecting_variant_count = intersection.count()
            _msg = f"Intersecting variants: {intersecting_variant_count}"
            print(_msg)
            self.log += _msg + '\n'
            
            # Filter to intersecting variants
            exome_filtered = exome_filtered.filter_rows(hl.is_defined(intersection[exome_filtered.row_key]))
            imputation_filtered = imputation_filtered.filter_rows(hl.is_defined(intersection[imputation_filtered.row_key]))
            
            total_intersecting_genotypes = intersecting_variant_count * n_samples
            _msg = f"Total intersecting genotypes: {total_intersecting_genotypes}"
            print(_msg)
            self.log += _msg + '\n'
            
            # Annotate with AC and MAF bins
            exome_filtered = exome_filtered.annotate_rows(
                ac_bin=hl.case()
                    .when(exome_filtered.AC <= 1, "1")
                    .when(exome_filtered.AC <= 2, "2")
                    .when(exome_filtered.AC <= 3, "3")
                    .when(exome_filtered.AC <= 4, "4")
                    .when(exome_filtered.AC <= 5, "5")
                    .when(exome_filtered.AC <= 10, "6-10")
                    .default("10+"),
                maf_bin=hl.case()
                    .when(exome_filtered.MAF <= 0.02, "1-2%")
                    .when(exome_filtered.MAF <= 0.05, "2-5%")
                    .default("5%+")
            )
            
            imputation_filtered = imputation_filtered.annotate_rows(
                ac_bin=hl.case()
                    .when(imputation_filtered.AC <= 1, "1")
                    .when(imputation_filtered.AC <= 2, "2")
                    .when(imputation_filtered.AC <= 3, "3")
                    .when(imputation_filtered.AC <= 4, "4")
                    .when(imputation_filtered.AC <= 5, "5")
                    .when(imputation_filtered.AC <= 10, "6-10")
                    .default("10+"),
                maf_bin=hl.case()
                    .when(imputation_filtered.MAF <= 0.02, "1-2%")
                    .when(imputation_filtered.MAF <= 0.05, "2-5%")
                    .default("5%+")
            )
        else:
            exome_filtered = self.exo
            imputation_filtered = self.imp
            n_samples = exome_filtered.count_cols()
            intersecting_variant_count = exome_filtered.count_rows()
            total_intersecting_genotypes = intersecting_variant_count * n_samples
        
        _msg = "Computing global concordance"
        print(_msg)
        self.log += _msg + '\n'
        
        global_conc, cols_conc, rows_conc = hl.concordance(exome_filtered, imputation_filtered)
        
        variant_concordance = rows_conc.select(rows_conc.concordance)
        
        # Get AC and MAF bins
        try:
            exome_bins = exome_filtered.rows().select('ac_bin', 'maf_bin')
            imputation_bins = imputation_filtered.rows().select('ac_bin', 'maf_bin')
        except Exception as e:
            _msg = f"Error getting AC and MAF bins: {str(e)}\nSkipping bin analysis."
            print(_msg)
            self.log += _msg + '\n'
            exome_bins = None
            imputation_bins = None
        
        ac_bins = ["1", "2", "3", "4", "5", "6-10", "10+"]
        maf_bins = ["1-2%", "2-5%", "5%+"]
        
        # Initialize result lists
        ac_results = []
        ac_results_imputation = []
        maf_results = []
        maf_results_imputation = []
        
        # Only perform bin analysis if bins are available
        if exome_bins is not None and imputation_bins is not None:
            # Process AC bins for exome
            for ac_bin in ac_bins:
                _msg = f"Processing AC bin {ac_bin} (exome)"
                print(_msg)
                self.log += _msg + '\n'
                
                bin_variants = exome_bins.filter(exome_bins.ac_bin == ac_bin)
                n_variants = bin_variants.count()
                
                if n_variants == 0:
                    _msg = f"No variants in AC bin {ac_bin} (exome)"
                    print(_msg)
                    self.log += _msg + '\n'
                    empty_conc = [[0 for _ in range(5)] for _ in range(5)]
                    ac_results.append({
                        'bin': ac_bin,
                        'concordance_matrix': empty_conc,
                        'n_variants': 0,
                        'het_genotypes_count': 0,
                        **self.compute_concordance_metrics(empty_conc)
                    })
                    continue
                
                het_count = exome_filtered.filter_rows(exome_filtered.ac_bin == ac_bin).aggregate_entries(
                    hl.agg.count_where(exome_filtered.GT.is_het())
                )
                
                bin_concordance = variant_concordance.semi_join(bin_variants)
                
                bin_conc_matrix = [[0 for _ in range(5)] for _ in range(5)]
                
                for i in range(5):
                    for j in range(5):
                        bin_conc_matrix[i][j] = bin_concordance.aggregate(
                            hl.agg.sum(bin_concordance.concordance[i][j])
                        )
                
                metrics = self.compute_concordance_metrics(bin_conc_matrix)
                
                ac_results.append({
                    'bin': ac_bin,
                    'concordance_matrix': bin_conc_matrix,
                    'n_variants': n_variants,
                    'het_genotypes_count': het_count,
                    **metrics
                })
            
            # Process AC bins for imputation
            for ac_bin in ac_bins:
                _msg = f"Processing AC bin {ac_bin} (imputation)"
                print(_msg)
                self.log += _msg + '\n'
                
                bin_variants = imputation_bins.filter(imputation_bins.ac_bin == ac_bin)
                n_variants = bin_variants.count()
                
                if n_variants == 0:
                    _msg = f"No variants in AC bin {ac_bin} (imputation)"
                    print(_msg)
                    self.log += _msg + '\n'
                    empty_conc = [[0 for _ in range(5)] for _ in range(5)]
                    ac_results_imputation.append({
                        'bin': ac_bin,
                        'concordance_matrix': empty_conc,
                        'n_variants': 0,
                        'het_genotypes_count': 0,
                        **self.compute_concordance_metrics(empty_conc)
                    })
                    continue
                
                het_count = imputation_filtered.filter_rows(imputation_filtered.ac_bin == ac_bin).aggregate_entries(
                    hl.agg.count_where(imputation_filtered.GT.is_het())
                )
                
                bin_concordance = variant_concordance.semi_join(bin_variants)
                
                bin_conc_matrix = [[0 for _ in range(5)] for _ in range(5)]
                
                for i in range(5):
                    for j in range(5):
                        bin_conc_matrix[i][j] = bin_concordance.aggregate(
                            hl.agg.sum(bin_concordance.concordance[i][j])
                        )
                
                metrics = self.compute_concordance_metrics(bin_conc_matrix)
                
                ac_results_imputation.append({
                    'bin': ac_bin,
                    'concordance_matrix': bin_conc_matrix,
                    'n_variants': n_variants,
                    'het_genotypes_count': het_count,
                    **metrics
                })
            
            # Process MAF bins for exome
            for maf_bin in maf_bins:
                _msg = f"Processing MAF bin {maf_bin} (exome)"
                print(_msg)
                self.log += _msg + '\n'
                
                bin_variants = exome_bins.filter(exome_bins.maf_bin == maf_bin)
                n_variants = bin_variants.count()
                
                if n_variants == 0:
                    _msg = f"No variants in MAF bin {maf_bin} (exome)"
                    print(_msg)
                    self.log += _msg + '\n'
                    empty_conc = [[0 for _ in range(5)] for _ in range(5)]
                    maf_results.append({
                        'bin': maf_bin,
                        'concordance_matrix': empty_conc,
                        'n_variants': 0,
                        'het_genotypes_count': 0,
                        **self.compute_concordance_metrics(empty_conc)
                    })
                    continue
                
                het_count = exome_filtered.filter_rows(exome_filtered.maf_bin == maf_bin).aggregate_entries(
                    hl.agg.count_where(exome_filtered.GT.is_het())
                )
                
                bin_concordance = variant_concordance.semi_join(bin_variants)
                
                bin_conc_matrix = [[0 for _ in range(5)] for _ in range(5)]
                
                for i in range(5):
                    for j in range(5):
                        bin_conc_matrix[i][j] = bin_concordance.aggregate(
                            hl.agg.sum(bin_concordance.concordance[i][j])
                        )
                
                metrics = self.compute_concordance_metrics(bin_conc_matrix)
                
                maf_results.append({
                    'bin': maf_bin,
                    'concordance_matrix': bin_conc_matrix,
                    'n_variants': n_variants,
                    'het_genotypes_count': het_count,
                    **metrics
                })
            
            # Process MAF bins for imputation
            for maf_bin in maf_bins:
                _msg = f"Processing MAF bin {maf_bin} (imputation)"
                print(_msg)
                self.log += _msg + '\n'
                
                bin_variants = imputation_bins.filter(imputation_bins.maf_bin == maf_bin)
                n_variants = bin_variants.count()
                
                if n_variants == 0:
                    _msg = f"No variants in MAF bin {maf_bin} (imputation)"
                    print(_msg)
                    self.log += _msg + '\n'
                    empty_conc = [[0 for _ in range(5)] for _ in range(5)]
                    maf_results_imputation.append({
                        'bin': maf_bin,
                        'concordance_matrix': empty_conc,
                        'n_variants': 0,
                        'het_genotypes_count': 0,
                        **self.compute_concordance_metrics(empty_conc)
                    })
                    continue
                
                het_count = imputation_filtered.filter_rows(imputation_filtered.maf_bin == maf_bin).aggregate_entries(
                    hl.agg.count_where(imputation_filtered.GT.is_het())
                )
                bin_concordance = variant_concordance.semi_join(bin_variants)
                bin_conc_matrix = [[0 for _ in range(5)] for _ in range(5)]
                
                for i in range(5):
                    for j in range(5):
                        bin_conc_matrix[i][j] = bin_concordance.aggregate(
                            hl.agg.sum(bin_concordance.concordance[i][j])
                        )
                
                metrics = self.compute_concordance_metrics(bin_conc_matrix)
                
                maf_results_imputation.append({
                    'bin': maf_bin,
                    'concordance_matrix': bin_conc_matrix,
                    'n_variants': n_variants,
                    'het_genotypes_count': het_count,
                    **metrics
                })
        
        # Compile results
        results = {
            'ac_results': ac_results,
            'ac_results_imputation': ac_results_imputation,
            'maf_results': maf_results,
            'maf_results_imputation': maf_results_imputation,
            'global_concordance': global_conc,
            'n_samples': n_samples,
            'intersecting_variant_count': intersecting_variant_count,
            'total_intersecting_genotypes': total_intersecting_genotypes
        }
        
        # Output results to JSON if requested
        if output_json:
            try:
                json_out_path = f'{self.out_path_prefix}concordance_results.json'
                with hl.hadoop_open(json_out_path, 'w') as f:
                    json.dump(results, f)
                _msg = f"Saved concordance results to {json_out_path}"
                print(_msg)
                self.log += _msg + '\n'
            except Exception as e:
                _msg = f"Error saving JSON results: {str(e)}"
                print(_msg)
                self.log += _msg + '\n'
        
        # Export TSV files if output_name is provided (backward compatibility)
        if output_name:
            try:
                if isinstance(output_name, list):
                    col_conc_fname = f'{self.out_path_prefix}{output_name[0]}'
                    cols_conc.export(col_conc_fname)
                    row_conc_fname = f'{self.out_path_prefix}{output_name[1]}'
                    rows_conc.export(row_conc_fname)
                elif isinstance(output_name, str):
                    col_conc_fname = f'{self.out_path_prefix}{output_name}_sample.tsv'
                    cols_conc.export(col_conc_fname)
                    row_conc_fname = f'{self.out_path_prefix}{output_name}_sites.tsv'
                    rows_conc.export(row_conc_fname)
                else:
                    raise ValueError("output_name must be either string or list.")
                    
                _msg = f"Saved sample and site concordance tables to {col_conc_fname} and {row_conc_fname}."
                print(_msg)
                self.log += _msg + '\n'
                
                # Store for backward compatibility
                self.gc = global_conc
                self.sample_conc = pd.read_csv(col_conc_fname, sep='\t', header=0)
                self.site_conc = pd.read_csv(row_conc_fname, sep='\t', header=0)
                
                # Generate histograms if requested
                if graph_out:
                    try:
                        self.graph_histogram(self.sample_conc['n_discordant'], 'sample')
                        self.graph_histogram(self.site_conc['n_discordant'], 'variant')
                    except Exception as e:
                        _msg = f"Error generating histograms: {str(e)}"
                        print(_msg)
                        self.log += _msg + '\n'
            except Exception as e:
                _msg = f"Error exporting concordance tables: {str(e)}"
                print(_msg)
                self.log += _msg + '\n'
        
        # Update the class with results
        self.concordance_results = results
        
        return results

    def graph_histogram(self, data, discordance_type):
        """
        Generate a histogram for concordance data.
        
        Args:
            data (array-like): Data to plot.
            discordance_type (str): Type of discordance ('sample' or 'variant').
        """
        if discordance_type not in {'sample', 'variant'}:
            raise ValueError("Please provide either sample or variant discordance.")
        
        try:
            hist, edges = np.histogram(data)
            p = figure(title=f"Distribution of {discordance_type} discordance", 
                      x_axis_label="n_discordant", 
                      y_axis_label="Frequency")
            p.quad(top=hist, bottom=0, left=edges[:-1], right=edges[1:], fill_color="navy", line_color="white", alpha=0.5)
            show(p)
        except Exception as e:
            _msg = f"Error generating histogram: {str(e)}"
            print(_msg)
            self.log += _msg + '\n'

    def create_summary_dataframes(self, results=None):
        """
        Create summary dataframes from concordance results.
        
        Args:
            results (dict, optional): Concordance results. If None, uses self.concordance_results.
            
        Returns:
            tuple: A tuple containing (ac_df, maf_df) - dataframes for AC and MAF results.
        """
        try:
            if results is None:
                if hasattr(self, 'concordance_results'):
                    results = self.concordance_results
                else:
                    raise ValueError("No concordance results available. Run concordance() first.")
            
            def process_results(results_list, source):
                df = pd.DataFrame([{
                    'bin': r['bin'],
                    'concordance_rate': r['concordance_rate'],
                    'nonref_concordance': r['nonref_concordance'],
                    'n_variants': r['n_variants'],
                    'n_discordant': r['n_discordant'],
                    'n_total': r['n_total'],
                    'het_to_ref': r['het_to_ref'],
                    'het_to_het': r['het_to_het'],
                    'het_genotypes_count': r['het_genotypes_count'],
                    'source': source
                } for r in results_list])
                return df
            
            ac_df_exome = process_results(results['ac_results'], 'exome')
            ac_df_imputation = process_results(results['ac_results_imputation'], 'imputation')
            ac_df = pd.concat([ac_df_exome, ac_df_imputation])
            
            ac_df['bin'] = pd.Categorical(ac_df['bin'], categories=['1', '2', '3', '4', '5', '6-10', '10+'], ordered=True)
            ac_df = ac_df.sort_values(['source', 'bin'])
            
            maf_df_exome = process_results(results['maf_results'], 'exome')
            maf_df_imputation = process_results(results['maf_results_imputation'], 'imputation')
            maf_df = pd.concat([maf_df_exome, maf_df_imputation])
            
            maf_df['bin'] = pd.Categorical(maf_df['bin'], categories=['1-2%', '2-5%', '5%+'], ordered=True)
            maf_df = maf_df.sort_values(['source', 'bin'])
            
            return ac_df, maf_df
        except Exception as e:
            import traceback
            _msg = f"Error creating summary dataframes: {str(e)}\n{traceback.format_exc()}"
            print(_msg)
            self.log += _msg + '\n'
            return None, None

    def display_result_tables(self, ac_df=None, maf_df=None, results=None):
        """
        Display formatted concordance result tables.
        
        Args:
            ac_df (DataFrame, optional): DataFrame with AC results.
            maf_df (DataFrame, optional): DataFrame with MAF results.
            results (dict, optional): Concordance results. Used if ac_df or maf_df is None.
            
        Returns:
            tuple: A tuple containing (ac_table, maf_table) - formatted tables.
        """
        try:
            if ac_df is None or maf_df is None:
                if results is None:
                    if hasattr(self, 'concordance_results'):
                        results = self.concordance_results
                    else:
                        raise ValueError("No concordance results available. Run concordance() first.")
                ac_df, maf_df = self.create_summary_dataframes(results)
            
            def format_table(df, bin_name):
                table = df.copy()
                table['concordance_rate'] = table['concordance_rate'].map('{:.3%}'.format)
                table['nonref_concordance'] = table['nonref_concordance'].map('{:.3%}'.format)
                table = table.rename(columns={
                    'bin': f'{bin_name} Bin',
                    'n_variants': '# Intersecting Variants',
                    'n_discordant': '# Discordant',
                    'n_total': '# Total Genotypes',
                    'concordance_rate': 'Concordance Rate',
                    'nonref_concordance': 'Non-ref Concordance',
                    'het_to_ref': 'Het→Ref',
                    'het_to_het': 'Het→Het',
                    'het_genotypes_count': 'Het Genotypes',  
                    'source': 'Source'
                })
                return table
            
            ac_table = format_table(ac_df, 'AC')
            _msg = "Allele Count (AC) Results:"
            print(_msg)
            self.log += _msg + '\n'
            
            for source in ['exome', 'imputation']:
                _msg = f"\n{source.upper()} AC Results:"
                print(_msg)
                self.log += _msg + '\n'
                
                source_table = ac_table[ac_table['Source'] == source].drop('Source', axis=1)
                _msg = source_table.to_string(index=False)
                print(_msg)
                self.log += _msg + '\n'
            
            maf_table = format_table(maf_df, 'MAF')
            _msg = "\nMinor Allele Frequency (MAF) Results:"
            print(_msg)
            self.log += _msg + '\n'
            
            for source in ['exome', 'imputation']:
                _msg = f"\n{source.upper()} MAF Results:"
                print(_msg)
                self.log += _msg + '\n'
                
                source_table = maf_table[maf_table['Source'] == source].drop('Source', axis=1)
                _msg = source_table.to_string(index=False)
                print(_msg)
                self.log += _msg + '\n'
            
            _msg = "\nSummary Statistics:"
            print(_msg)
            self.log += _msg + '\n'
            
            _msg = f"Total Samples: {results['n_samples']}"
            print(_msg)
            self.log += _msg + '\n'
            
            _msg = f"Total Intersecting Variants: {results['intersecting_variant_count']}"
            print(_msg)
            self.log += _msg + '\n'
            
            _msg = f"Total Intersecting Genotypes: {results['total_intersecting_genotypes']}"
            print(_msg)
            self.log += _msg + '\n'
            
            _msg = f"Genotypes per variant: {results['n_samples']}"
            print(_msg)
            self.log += _msg + '\n'
            
            return ac_table, maf_table
        except Exception as e:
            import traceback
            _msg = f"Error displaying result tables: {str(e)}\n{traceback.format_exc()}"
            print(_msg)
            self.log += _msg + '\n'
            return None, None

    def save_concordance_results(self, file_path=None):
        """
        Save concordance results to a JSON file.
        
        Args:
            file_path (str, optional): Path to save results. If None, uses default path.
            
        Returns:
            str: Path to the saved file, or None if there was an error.
        """
        try:
            if not hasattr(self, 'concordance_results'):
                raise ValueError("No concordance results available. Run concordance() first.")
            
            if file_path is None:
                file_path = f'{self.out_path_prefix}concordance_results.json'
            
            with hl.hadoop_open(file_path, 'w') as f:
                json.dump(self.concordance_results, f)
            
            _msg = f"Saved concordance results to {file_path}"
            print(_msg)
            self.log += _msg + '\n'
            
            return file_path
        except Exception as e:
            import traceback
            _msg = f"Error saving concordance results: {str(e)}\n{traceback.format_exc()}"
            print(_msg)
            self.log += _msg + '\n'
            return None

    def load_concordance_results(self, file_path=None):
        """
        Load concordance results from a JSON file.
        
        Args:
            file_path (str, optional): Path to load results from. If None, uses default path.
        
        Returns:
            dict: The loaded concordance results, or None if there was an error.
        """
        try:
            if file_path is None:
                file_path = f'{self.out_path_prefix}concordance_results.json'
            
            with hl.hadoop_open(file_path, 'r') as f:
                results = json.load(f)
            
            self.concordance_results = results
            
            _msg = f"Loaded concordance results from {file_path}"
            print(_msg)
            self.log += _msg + '\n'
            
            return results
        except Exception as e:
            import traceback
            _msg = f"Error loading concordance results: {str(e)}\n{traceback.format_exc()}"
            print(_msg)
            self.log += _msg + '\n'
            return None
            
    def remove_discordant(self, which, thresh=None):
        """
        Filter the merged dataset to remove samples or variants with high discordance.
        
        Args:
            which (str): Type of filtering ('sample', 'variant', or 'site').
            thresh (int): Threshold for number of discordant genotypes to allow.
        """
        if self.sample_conc is None or self.site_conc is None:
            raise ValueError('Please compute concordance tables first!')
        if not thresh:
            raise ValueError('Please provide a threshold.')
        if which not in set(['sample', 'variant', 'site']):
            raise ValueError('Must be either variant, site, or sample discordance.')
        
        if which == 'sample':
            df = self.sample_conc
        else:
            df = self.site_conc
        filtered = df[df['n_discordant'] > thresh]
        
        if which == 'site' or which == 'variant':
            set_keep = hl.literal(set(list(filtered['locus'].apply(lambda x: hl.parse_locus(x)))))
            _msg = f'Filtering {which} to < {thresh} n_discordant.'
            print(_msg)
            self.log += _msg + '\n'
            self.merged = self.merged.filter_rows(set_keep.contains(self.merged.locus), keep=True)
        else:
            set_keep = hl.literal(set(filtered['s']))
            _msg = f'Filtering {which} to < {thresh} n_discordant.'
            print(_msg)
            self.log += _msg + '\n'
            self.merged = self.merged.filter_cols(set_keep.contains(self.merged.s), keep=True)
        
    def calc_overlap(self, output_name=None):
        """
        calculates the overlap matrices of variants and samples between datasets, first filtering to bi-allelic sites
        Args:
            output_name (str, optional): Output file name.
        """
        self.overlap = self.exo.semi_join_rows(self.imp.rows())
        self.overlap = self.overlap.semi_join_cols(self.exo.cols())
        self._n_overlap_sites, self._n_overlap_samples = self.overlap.count()
        print(f'{self._n_overlap_sites} sites and {self._n_overlap_samples} samples are similar between {self.exo_path} and {self.imp_path}.')
        
        if output_name:
            out_name = f'{self.out_path_prefix}{output_name}.txt'
            locs = self.overlap.locus.collect()
            loci = list(map(lambda x: str(x), locs))
            self.export_flat(loci, out_name)
            print(f'{len(loci)} variants overlap between both files, written to {out_name}.')
            
    def plot_concordance(self, col_conc_fname):
        """
        Plot concordance vs MAF.
        
        Args:
            col_conc_fname (str): Path to column concordance file.
        """
        self.merged = hl.variant_qc(self.merged)
        self.merged = self.merged.annotate_rows(MAF=hl.min(self.merged.variant_qc.AF))
        
        col_conc = pd.read_csv(col_conc_fname, sep='\t', header=0)
        x = self.merged.MAF.collect()
        y = (np.ones(len(col_conc['n_discordant'].values)) * self.shape[1] - col_conc['n_discordant'].values) / self.shape[1]        
        
        z = np.array(list(sorted(zip(x, y), key=lambda v: v[0]))).T
        
        p = figure(title="Concordance over all chromosomes", 
                   x_axis_label='MAF (%)', 
                   y_axis_label='Sample concordance')
        self.x = x
        self.y = y
    
        p.circle(z[0, :]*100, z[1, :], size=3, alpha=0.5)
        show(p)
        self.maf_plot = p

    def generate_sites_file(self, export=False):
        """
        Generate a sites file with marker information.
        
        Args:
            export (bool, optional): If True, export the sites file.
            
        Returns:
            str: Generated sites string.
        """
        rsids = self.merged.rsid.collect()
        locus = self.merged.locus.collect()
        
        sites = 'MarkerName\tChromosome\tPosition\n'
        for rsid, loc in zip(rsids, locus):
            sites += f'{rsid}\t{loc.contig}\t{loc.position}\n'
            
        if export:
            export_outpath = f'{self.out_path_prefix}.sites'
            self.export_flat(sites, export_outpath)
            
        return sites
    
    def merge(self):
        """
        Merge exome and imputation datasets with priority.
        """
        def align_mt2_cols_to_mt1(mt1, mt2):
            mt1 = mt1.add_col_index()
            mt2 = mt2.add_col_index()
            new_col_order = mt2.index_cols(mt1.col_key).col_idx.collect()
            return mt2.choose_cols(new_col_order)

        if self.priority == 'exome':
            not_exo = self.imp.anti_join_rows(self.exo.rows())
            self.exo = self.exo.semi_join_cols(not_exo.cols())
            not_exo = align_mt2_cols_to_mt1(self.exo, not_exo)
            self.merged = self.exo.union_rows(not_exo)
        else:
            not_imp = self.exo.anti_join_rows(self.imp.rows())
            self.imp = self.imp.semi_join_cols(not_imp.cols())
            not_imp = align_mt2_cols_to_mt1(self.imp, not_imp)
            self.merged = self.imp.union_rows(not_imp)
    
    def save_log(self):
        """
        Saves log string to logging file.
        """
        self.export_flat(self.log, self.log_outpath)
    
    @property
    def shape(self):
        if self._shape is None:
            self._shape = self.merged.count()
        return self._shape
     
