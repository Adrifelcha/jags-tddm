load dic
load tddm
model in "test/tddm.bug"
data in "test/tddm.data"
compile, nchains(1)
initialize
update 5
monitor set drift, thin(1)
monitor set bound, thin(1)
monitor set theta0, thin(1)
monitor set ter0, thin(1)
monitor deviance
update 5
coda *, stem('samples_1_')

