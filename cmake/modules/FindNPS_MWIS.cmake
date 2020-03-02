# - Try to find NPS_MWIS
# Once done this will define
#  NPS_MWIS_FOUND - System has NPS_MWIS
#  NPS_MWIS_SRC_DIR
#  NPS_MWIS_INCLUDE_DIR

if (NPS_MWIS_SRC_DIR)
    # in cache already
    set(NPS_MWIS_FOUND TRUE)
else (NPS_MWIS_SRC_DIR)

    find_path(NPS_MWIS_INCLUDE_DIR
            NAMES algorithm.h
            PATHS ${CMAKE_CURRENT_SOURCE_DIR}/extern/nps_mwis/src)

    set(NPS_MWIS_SRC_DIR "${NPS_MWIS_INCLUDE_DIR}")

    include(FindPackageHandleStandardArgs)
    # handle the QUIETLY and REQUIRED arguments and set NPS_MWIS_FOUND to TRUE
    # if all listed variables are TRUE
    find_package_handle_standard_args(NPS_MWIS  DEFAULT_MSG
            NPS_MWIS_SRC_DIR)

    mark_as_advanced(NPS_MWIS_SRC_DIR)

endif(NPS_MWIS_SRC_DIR)
