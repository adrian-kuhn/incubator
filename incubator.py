#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
incubator - CLI application interface

Usage:
    incubator <setting-file>


Arguments:
    setting-file    The name of the setting file being loaded for the next run

Examples:
    incubator settings_watershed
"""

import sys
import logging.config
from docopt import docopt

from incubator import settings, processing, algorithms

log = logging.getLogger()


def load_assignments():
    """
    Create an assignment for every raw data dictionary object found in configured assignments settings.

    :return: List of assignments
    :rtype: List class:`Assignment`
    """
    assignments = []
    for assignment in settings.get("assignments"):
        assignments.append(processing.Assignment(processing.TreeTops.from_feature_class(assignment["target"]),
                                                 **assignment["raw_data"]))
    return assignments


def main():
    """
    Main function checks for python 3, initialize the settings and logging and starts the incubator process.
    """
    if sys.version_info.major != 3:
        print('Please activate the conda package and run on python 3')
        return

    args = docopt(__doc__, version='8.5.0')
    settings.init(args["<setting-file>"])
    logging.config.dictConfig(settings.get("logging"))

    try:
        algorithm = getattr(algorithms, settings.get("processing.algorithm"))
    except AttributeError:
        log.error("[] not found in module algorithms".format(settings.get("processing.algorithm")))
        return

    generations = settings.get("processing.number_of_generations")
    individuals = settings.get("processing.number_of_individuals")
    run_parallel = settings.get("processing.parallel")
    p = processing.Incubator(algorithm, load_assignments(), generations, individuals, run_parallel)
    p.breed()


if __name__ == "__main__":
    main()
