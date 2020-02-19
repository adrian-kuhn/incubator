#!/usr/bin/env python
# -*- coding: utf-8 -*-

import sys
import logging.config

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
        assignments.append(processing.Assignment(processing.TreePoints.from_feature_class(assignment["target"]),
                                                 **assignment["raw_data"]))
    return assignments


def main():
    """
    Main function checks for python 3, initialize the settings and logging and starts the incubator process.
    """
    if sys.version_info.major != 3:
        print('Please activate the conda package and run on python 3')
        return

    settings.init()
    logging.config.dictConfig(settings.get("logging"))

    # setup_las_files()

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
