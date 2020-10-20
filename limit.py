def get_limit( npzfile, CL=0.95, forcePositive=False ):
  import scipy.stats
  import numpy as np
  #convert 1-sided to 2-sided interval and calculate deltaTS value
  #See slide 33 in https://fermi.gsfc.nasa.gov/science/mtgs/summerschool/2013/week1/ML_intro.pdf
  deltaTS = scipy.stats.chi2.ppf(1.0-2*(1.0-CL),1)
  try:
  #if True:
    res = np.load(npzfile)
    a = res["a"]
    ll = res["cc"]
    a = a[ll<1e12]
    ll = ll[ll<1e12]
    #find minimum
    imin = np.argmin( ll )
    if a[imin] < 0 and forcePositive:
      imin = np.argmin(np.abs(a))
    target = ll[imin] + deltaTS/2.0
    print imin, ll[imin], target
    for i in range(0, imin):
      if ( ll[i] <= target and ll[i+1] >= target ) or ( ll[i] >= target and ll[i+1] <= target ):
        m = (ll[i+1] - ll[i]) / (a[i+1]-a[i])
        val1 = (a[i] + (target - ll[i]) / m )
        print "Target1 (%.1f) crossed between %.2f and %.2f; linear interpolation: %.2f." %( target, a[i], a[i+1], val1)
    return val1
  except:
    return np.nan  
