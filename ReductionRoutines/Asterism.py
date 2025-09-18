from astrom import *
from astrom import \
     find_asterism_star_locations, \
     comparison
import configparser

# configuration data
config = configparser.ConfigParser() # for parsing values in .init file
config.read("astrom/config.ini")

dateStringShort = config["dataset_string"]["DATE_SHORT"]

# make the directories
make_dirs()


### PART 1: FIND DEWARP SOLUTION
# match the empirical and ideal pinholes
# barrelCenterPass is (x,y)
find_asterism_star_locations.find_stars(
    dateString = dateStringShort,
    number_of_dithers = 17)


# find position angle offset and plate scale
comparison.angOffset_plateScale(dateStringShort,
                                config["dataset_string"]["ASTERISM_PLOT_TITLE_STRING"],
                                plot=True)