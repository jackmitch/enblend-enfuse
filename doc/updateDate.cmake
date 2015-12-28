file(STRINGS ${INPUT_PATH}/${DATEFILE} UPDATE_DATE)
configure_file(${INPUT_PATH}/version.in "${INPUT_PATH}/${FILENAME}")
