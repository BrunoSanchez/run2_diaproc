from lsst.dia.pipe.selectImages import TimeBestSeeingWcsSelectImagesTask
print(config)
config.select.retarget(TimeBestSeeingWcsSelectImagesTask)

config.select.minMJD = 59215.0 # jan 1 2021
#config.select.minMJD = 59580.0  # jan 1 2022
config.select.maxMJD = 59945.0  # jan 1 2023
#config.select.maxMJD = 60310.0  # jan 1 2024
#config.select.maxMJD = 61041.0  # jan 1 2026
#config.select.maxMJD = 61406.0  # jan 1 2026

config.select.minPsfFwhm = 2.5
config.select.maxPsfFwhm = 5.0


