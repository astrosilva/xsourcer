from target import Target

source = Target('CXOGC_J174543.4-285841')
source.query()
source.download(242, 3392, 3665, 13438, 22707, mode='q')
#source.reprocess()
source.merge_obs()
source.periodogram()