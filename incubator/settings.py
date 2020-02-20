#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Module for loading the settings for a given environment from json file
"""
import json
import os

__SETTINGS = None
""" Global settings """

__SETTING_DIR = "settings"
""" The folder for setting files """


def init(settings_file):
    """
    Search the settings file in settings directory and given environment variable.
    The settings file has to be in the structure 'env'.settings.json
    Initialize the settings by loading the JSON file and return the settings Dictionary.

    :param settings_file: Name of the setting file. Default is 'settings'
    :type settings_file: String
    :raises SettingsNotFoundException: Raises a SettingsNotFoundException, if the settings file cannot be found.
    """
    global __SETTINGS

    try:
        settings_file = "{}.json".format(settings_file) if not settings_file.endswith(".json") else settings_file
        s = open(os.path.join(__SETTING_DIR, settings_file), 'rt')
        __SETTINGS = json.load(s)

    except (OSError, IOError, FileNotFoundError):
        raise SettingsNotFoundException(__SETTING_DIR, settings_file)


def get(key, settings=None):
    """
    Get a setting value by key. Settings can be nested and loaded by point (.) notation.

    :param key: The name of the setting variable
    :type key: String
    :param settings: The setting dictionary
    :type settings: Dictionary
    :return: The setting value.
    :rtype: String
    """
    settings = __SETTINGS if settings is None else settings
    parts = key.split(".", 1)
    if len(parts) == 2:
        return get(parts[1], get(parts[0], settings))

    if key not in settings:
        raise KeyError("Could not find {0} in settings. "
                       "Check initialization and settings variables".format(key))
    return settings[key]


class SettingsNotFoundException(Exception):
    """
    Custom exception for invalid settings file. Includes the provided environments in the exception message
    """
    def __init__(self, settings_dir, setting_file_postfix):
        """
        Constructor needs the base path to the setting files and the required setting file postfix.

        :param settings_dir: Absolute path to the settings folder
        :type settings_dir: String
        :param setting_file_postfix: Required postfix
        :type settings_dir: String
        """
        environments = [j.replace(setting_file_postfix, "") for j in os.listdir(settings_dir)]
        message = "Setting file not found in {0}. " \
                  "Call settings.init(settings_path, env) to use tree detection in the preferred environment. " \
                  "Found environment settings for {1}".format(settings_dir, ', '.join(environments))
        super().__init__(message)
