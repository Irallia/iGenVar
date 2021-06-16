from argparse import ArgumentParser
import csv
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt


class Callset(object):
    def __init__(self, caller, precision, recall, opacity):
        self.caller = caller
        self.precision = float(precision)
        self.recall = float(recall)
        self.opacity = float(opacity)
        self.shape = None
        self.colour = None

    def plot(self):
        return plt.scatter(x=self.recall,
                           y=self.precision,
                           c=self.colour,
                           marker=self.shape,
                           alpha=self.opacity)


def main():
    args = get_args()
    calls = [Callset(*line) for line in csv.reader(open(args.results), delimiter="\t")]
    caller_index, callers = caller_to_colour(calls)
    lines = [c.plot() for c in calls]
    plt.xlim(0, args.max)
    plt.ylim(0, args.max)
    plt.xlabel('Recall')
    plt.ylabel('Precision')
    if len(callers) > 1:
        plt.gca().add_artist(plt.legend(handles=[lines[i] for i in caller_index],
                                        labels=callers,
                                        loc='lower {}'.format(args.legend),
                                        frameon=False))
    plt.tight_layout()
    plt.savefig("precision-recall.png")


def caller_to_colour(calls):
    """
    Assign colours to different callers in the input file
    Set the attribute of the class instances

    return a list of indices for which each caller is found uniquely and all callers
    sorted by callers
    """
    colours = plt.rcParams['axes.prop_cycle'].by_key()['color'] * 10
    callers = sorted(set([c.caller for c in calls]), reverse=True)
    caller_to_colour_dict = {ca: co for ca, co in zip(callers, colours)}
    for c in calls:
        c.colour = caller_to_colour_dict[c.caller]
    index_and_callers = zip([[c.caller for c in calls].index(i) for i in callers], callers)
    return zip(*sorted(index_and_callers, key=lambda x: x[1]))


def get_args():
    parser = ArgumentParser(description="Plot precision and recall for multiple callsets")
    parser.add_argument("results",
                        help="A tsv file with caller precision recall")
    parser.add_argument("--max",
                        help="A maximum value to which precision recall plot has to be limited",
                        type=int,
                        default=100)
    parser.add_argument("--legend",
                        help="Put legend left or right",
                        default="left")
    return parser.parse_args()


if __name__ == '__main__':
    main()
