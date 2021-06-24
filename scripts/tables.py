import os


def save_table(df, output_dir, filename):
    results_dir = f'tables/{output_dir}'
    if not os.path.isdir(results_dir):
        os.makedirs(results_dir)
    df.to_csv(results_dir + '/' + filename, float_format='%.3f')


def generate_tables(df, output_dir):
    df = df.reset_index()
    df = df.set_index(['date', 'region'])
    df = df.sum(level='region')
    df['proportion_sequenced'] = df.sequences / df.cases
    save_table(df, output_dir, 'region_table.csv')
