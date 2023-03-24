SUBDIRS = distributions

jagsmod_LTLIBRARIES = tddm.la

tddm_la_SOURCES = tddm.cc

tddm_la_CPPFLAGS = -I$(includedir)/JAGS

tddm_la_LIBADD = distributions/tddmdist.la
tddm_la_LIBADD += -ljags -ljrmath

tddm_la_LDFLAGS = -module -avoid-version



