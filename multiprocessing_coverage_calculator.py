import subprocess
import pandas as pd
import io
import os
from tqdm import tqdm
import multiprocessing


def process_gene(args):
    bam_file, ref, gene = args
    chrom = gene['Chromosome']
    start = gene['Start']
    end = gene['End']
    gene_id = gene['gene_id']
    strand = gene['Strand']

    pysamstats_command = f"pysamstats --type coverage_strand --min-baseq 0 --no-dup --no-del -c {chrom} -s {start} -e {end} {bam_file}"
    output = subprocess.check_output(pysamstats_command, shell=True, text=True)
    df_coverage = pd.read_csv(io.StringIO(output), sep='\t')
    df_coverage = df_coverage[(df_coverage['pos'] >= start) & (df_coverage['pos'] <= end)]

    return {
        'Gene_ID': gene_id,
        'Strand': strand,
        'all_reads': list(df_coverage['reads_all']),
        'forward_reads': list(df_coverage['reads_fwd']),
        'reverse_reads': list(df_coverage['reads_rev'])
    }

def calculate_gene_coverage_parallel(bam_file, ref):
    subprocess.run(["samtools", "index", bam_file])
    coverage_data = []

    # Create a pool of worker processes
    pool = multiprocessing.Pool()

    # Prepare arguments for each gene
    gene_args = [(bam_file, ref, gene) for _, gene in ref.iterrows()]

    # Process genes in parallel
    results = list(tqdm(pool.imap(process_gene, gene_args), total=len(gene_args)))

    # Close the pool to free up resources
    pool.close()
    pool.join()

    coverage_data.extend(results)

    df = pd.DataFrame(coverage_data)
    #df['Strand_specific_depth'] = df.apply(lambda row: row['reverse_reads'] if row['Strand'] == '-' else row['forward_reads'], axis=1)
    df = df[['Gene_ID', 'all_reads']]
    

    return df
    
    
if __name__ == '__main__':
    input_directory = '/lustre/abuzar.khan/Coverage_project/bamfile/pysamstat_cov/newbamfile'
    gtf_file_path = os.path.join('/lustre/abuzar.khan/Coverage_project/bamfile/pysamstat_cov/newbamfile', 'selected_gtf_for_51_genes.csv')
    df = pd.read_csv(gtf_file_path,index_col = 0)
    bam_files = [os.path.join(input_directory, file_name) for file_name in os.listdir(input_directory) if file_name.endswith(".bam")]
    #bam_files =bam_files[0:2]
    save_directory = '/lustre/abuzar.khan/Coverage_project/bamfile/pysamstat_cov/newbamfile'
    for bamfile in bam_files:
        output_csv_file = os.path.join(save_directory, 'Coverage', os.path.basename(bamfile).split('.')[0] + 'plot_dynamics_51_genes_non_specific.csv')
        cov_df = calculate_gene_coverage_parallel(bamfile,df)
        cov_df.to_csv(output_csv_file, index=False)
        


