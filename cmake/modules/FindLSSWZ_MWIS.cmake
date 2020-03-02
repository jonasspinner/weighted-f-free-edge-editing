# - Try to find LSSWZ_MWIS
# Only try to find the weighted_ls executable.
# Once done this will define
#  LSSWZ_MWIS_FOUND - System has LSSWZ_MWIS
#  LSSWZ_MWIS_SCRIPT - The LSSWZ_MWIS weighted_ls path

if (LSSWZ_MWIS_SCRIPT)
    # in cache already
    set(LSSWZ_MWIS_FOUND TRUE)
else (LSSWZ_MWIS_SCRIPT)

    find_file(LSSWZ_MWIS_SCRIPT
            NAMES weighted_ls
            PATHS ${CMAKE_CURRENT_SOURCE_DIR}/extern/lsswz_mwis/code/build)

    include(FindPackageHandleStandardArgs)
    # handle the QUIETLY and REQUIRED arguments and set LSSWZ_MWIS_FOUND to TRUE
    # if all listed variables are TRUE
    find_package_handle_standard_args(LSSWZ_MWIS  DEFAULT_MSG
            LSSWZ_MWIS_SCRIPT)

    mark_as_advanced(LSSWZ_MWIS_SCRIPT)

endif(LSSWZ_MWIS_SCRIPT)
