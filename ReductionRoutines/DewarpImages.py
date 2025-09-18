from astrom import *
from astrom import \
     match_pinholes, \
     find_dewarp_solution, \
     apply_dewarp_solution
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
match_pinholes.match_pinholes(
    translationPass = [0,-10],
    holeSpacingPass = 48.0,
    barrelCenterPass = [1124,715],
    barrelAmountPass = 0,
    rotationAnglePass = -0.8,
    writeoutString = config["dataset_string"]["DATASET_STRING"],
    plotTitleString = config["dataset_string"]["PINHOLE_PLOT_TITLE_STRING"],
    plot=True) # add another flag for when user is just testing model grid first


# work out the solution
find_dewarp_solution.find_dewarp(
    fileString = config["dataset_string"]["DATASET_STRING"],
    dateString = config["dataset_string"]["PINHOLE_PLOT_TITLE_STRING"],
    plot=True)


### PART 2: APPLY DEWARP
# apply dewarp solution to asterism frames
apply_dewarp_solution.apply_dewarp(
    writeoutString = config["dataset_string"]["DATASET_STRING"],
    maskUnsampled = True)