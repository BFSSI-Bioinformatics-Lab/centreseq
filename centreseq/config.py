import shutil
import logging

main_log = logging.getLogger('main_log')

EXTERNAL_DEPENDENCIES = ['mmseqs', 'prokka', 'muscle', 'seqkit', 'snp-sites']


def dependency_check(dependency: str) -> bool:
    """
    Checks if a given program is present in the user's $PATH
    :param dependency: String of program name
    :return: True if program is in $PATH, False if not
    """
    check = shutil.which(dependency)
    if check is not None:
        return True
    else:
        return False


def check_dependencies() -> bool:
    # Dependency check
    main_log.info("Conducting dependency check...")
    dependency_dict = dict()
    for dependency in EXTERNAL_DEPENDENCIES:
        dependency_dict[dependency] = dependency_check(dependency)
    if False in dependency_dict.values():
        main_log.error("ERROR: Cannot locate some dependencies in $PATH...")
        for key, value in dependency_dict.items():
            if not value:
                main_log.error(f"Dependency missing: {key}")
        return False
    else:
        for key, value in dependency_dict.items():
            main_log.debug(f"Dependency {key}: {value}")
    main_log.info("Dependencies OK")
    return True
