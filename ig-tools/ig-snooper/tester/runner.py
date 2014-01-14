from itertools import product, zip_longest
from IPython.parallel import Client, require
from utils import log_result


def main():
    client = Client()
    # view = client.load_balanced_view()
    view = client[:]

    @view.parallel()
    @require('os')
    @require('shutil')
    @require('utils')
    @require('itertools')
    def run_experiment(parameter_sets):
        REPEAT_COUNT = 10
        results = []
        for ml_window, avg_window, merge_threshold in parameter_sets:

            dir = utils.TEMP_DIR + str(os.getpid())
            if not os.path.exists(dir):
                os.mkdir(dir)
            cumulative_errs = ()
            for _ in range(REPEAT_COUNT):
                errs = utils.experiment(os.path.join(utils.DATA_DIR, 'vh.fasta'), os.path.join(utils.DATA_DIR, 'vh.kabat'), ml_window, avg_window, merge_threshold, 2000, dir)
                cumulative_errs = tuple(sum(es) for es in itertools.zip_longest(cumulative_errs, errs, fillvalue=0))
            # shutil.rmtree(dir)

            # have to use loop instead of list comprehension, because embedding REPEAT_COUNT in a comprehension results
            # in 'cannot pickle code object with closures' error
            final_errs = []
            for e in cumulative_errs:
                final_errs.append(e / REPEAT_COUNT)
            results.append(tuple(final_errs))

            with open(os.path.join(dir, 'log.txt'), 'a+') as log:
                utils.log_result(log, ml_window, avg_window, merge_threshold, *final_errs)

        return results


    parameter_sets = list(product(range(3, 34)[::2], range(3, 34)[::2], range(2, 10)))
    # parameter_sets = list(product([11, 13, 15], [3], [4, 5, 6, 7, 8, 9]))
    # parameter_sets = [(5, 5, 5)]
    # with open('missing.txt') as missing:
        # parameter_sets = [(tuple(n for n in map(float, line.split()))) for line in missing]
    async = run_experiment(parameter_sets)

    with open('log.txt', 'w') as log:
        for parameters, errs in zip_longest(parameter_sets, async):
            ml_window, avg_window, merge_threshold = parameters
            log_result(log, ml_window, avg_window, merge_threshold, *errs)


if __name__ == '__main__':
    main()