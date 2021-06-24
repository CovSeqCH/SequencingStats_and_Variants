#! /usr/bin/env python3

import argparse
import datetime


def valid_date(s):
    try:
        return datetime.datetime.strptime(s, "%Y-%m-%d")
    except ValueError:
        msg = "Not a valid date: '{0}'.".format(s)
        raise argparse.ArgumentTypeError(msg)


def main(args):
    """ Main entry point of the app """

    if(args.download):
        import scripts.download_sequences
        scripts.download_sequences.generate_csv()
        print('data updated')

    if(args.plot):
        import scripts.plots
        df = scripts.plots.load_data()
        df = scripts.plots.restrict_dates(df, args.start_date, args.end_date)
        scripts.plots.generate_plots(df, args.output_dir)
        print('plots generated')

    if(args.table):
        import scripts.plots
        import scripts.tables
        df = scripts.plots.load_data()
        df = scripts.plots.restrict_dates(df, args.start_date, args.end_date)
        scripts.tables.generate_tables(df, args.output_dir)
        print('tables generated')

    if(not any([args.download, args.plot, args.table])):
        parser.print_help()


if __name__ == "__main__":
    """ This is executed when run from the command line """
    parser = argparse.ArgumentParser()

    parser.add_argument("-d", "--download",
                        action="store_true", help="download data, if unset, use stored data")

    parser.add_argument("-p", "--plot",
                        action="store_true", help="create plots")

    parser.add_argument("-t", "--table",
                        action="store_true", help="create tables")

    parser.add_argument("-s",
                        "--start-date",
                        help="start date - format YYYY-MM-DD",
                        required=False,
                        type=valid_date)

    parser.add_argument("-e",
                        "--end-date",
                        help="end date - format YYYY-MM-DD",
                        required=False,
                        type=valid_date)

    parser.add_argument("-o",
                        "--output-dir",
                        help="output sub-directory",
                        required=False,
                        default='default')

    args = parser.parse_args()
    main(args)
