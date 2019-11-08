from lsst.dia.pipe.selectImages import TimeBestSeeingWcsSelectImagesTask
print(config)
config.select.retarget(TimeBestSeeingWcsSelectImagesTask)

config.select.minMJD = 59215.0 # jan 1 2021
#config.select.minMJD = 59580.0  # jan 1 2022
config.select.maxMJD = 59945.0  # jan 1 2023
#config.select.maxMJD = 60310.0  # jan 1 2024
#config.select.maxMJD = 61041.0  # jan 1 2026
#config.select.maxMJD = 61406.0  # jan 1 2026

# scale plate should be 0.2 arcsecs/px
config.select.minPsfFwhm = 2.5  # = 0.5 arcsecs
config.select.maxPsfFwhm = 5.0  # = 1.0 arcsecs


