import sys

try:
    import numpy
    have_numpy=True
except:
    have_numpy=False

try:
    import scipy
    from scipy.interpolate import interp1d
    import scipy.weave
    have_scipy=True
except:
    have_scipy=False





def histogram(data, bins=None, binsize=1., min=None, max=None, rev=False, use_weave=True, verbose=0):
    """
    Similar to IDL histogram.

    For reverse indices, the fast version uses weave from scipy. This is the
    default.  If scipy is not available a slower version is used.

    """

    if not have_scipy:
        use_weave=False
     
    dowhere=False
    s = data.argsort()
    if min is not None:
        dmin = min
        dowhere=True
    else:
        dmin = data[s[0]]

    if max is not None:
        dmax = max
        dowhere=True
    else:
        dmax = data[s[-1]]

    bsize = float(binsize)
#    except TypeError :
#        print '[Warning] TypeError reached - type {0}'.format(type(binsize[0]))
#        bsize = binsize
#    finally : pass

    if dowhere:
        # where will preserve order, so subscript with s
        w,=numpy.where( (data[s] >= dmin) & (data[s] <= dmax) )
        if w.size == 0:
            raise ValueError("No data in specified min/max range\n")
        s = s[w]

    nbin = numpy.int64( (dmax-dmin)/bsize ) + 1

    if bins is not None :
        s = bins
        nbin = len(s)

    if rev:
        revsize = s.size + nbin+1
        revind = numpy.zeros(revsize, dtype='i8')
    else:
        # this is just a dummy variable
        revind=numpy.zeros(1, dtype='i8')
    hist = numpy.zeros(nbin, dtype='i8')

    # populate the array from nbin+1:nbin+1+s.size
    # with the sort indices.  Simultaneosly record bin
    # edges at the beginning of reverse indices

    if use_weave:
        _weave_dohist(data, s, bsize, hist, revind, dorev=rev,verbose=verbose)
    else:
        _dohist(data, s, bsize, hist, revind, dorev=rev)

    if rev:
        return hist, revind
    else:
        return hist


def _weave_dohist(data, s, binsize, hist, rev, dorev=False, verbose=0):

    if dorev:
        dorev=1
    else:
        dorev=0
    """
    Weave version of histogram with reverse_indices

    s is an index into data, sorted and possibly a subset
    """
    code = """

    int64_t nbin = hist.size();
    int64_t binnum_old = -1;

    // index of minimum value
    int64_t imin = s(0);
    for (int64_t i=0; i<s.size(); i++) {

        int64_t offset = i+nbin+1;
        int64_t data_index = s(i);


        if (dorev) {
            rev(offset) = data_index;
        }

        int64_t binnum = (int64_t) ( (data(data_index)-data(imin))/binsize);

        if (binnum >= 0 && binnum < nbin) {
            if (dorev && (binnum > binnum_old) ) {
                int64_t tbin = binnum_old + 1;
                while (tbin <= binnum) {
                    rev(tbin) = offset;
                    tbin++;
                }
            }
            hist(binnum) = hist(binnum) + 1;
            binnum_old = binnum;
        }
    }

    int64_t tbin = binnum_old + 1;
    while (tbin <= nbin) {
        if (dorev) {
            rev(tbin) = rev.size();
        }
        tbin++;
    }

    """

    scipy.weave.inline(code, ['data','s','binsize','hist','rev','dorev'],
                       type_converters = scipy.weave.converters.blitz, verbose=verbose)
    return


def _dohist(data, s, binsize, hist, revind, dorev=False):

    nbin=hist.size
    offset = nbin+1
    i=0
    binnum_old = -1

    dmin = data[s[0]]
    while i < s.size:
        data_index = s[i]
        if dorev:
            revind[offset] = data_index

        val = data[data_index]

        binnum = numpy.int64( (val-dmin)/binsize )
        #print 'binnum:',binnum,' binnum old:',binnum_old, 'val:',val
        if binnum >= 0 and binnum < nbin:
        #if binnum >= 0:
            if binnum > binnum_old:
                tbin = binnum_old + 1
                while tbin <= binnum:
                    if dorev:
                        revind[tbin] = offset
                        #print '\t\trevind[%d] = %d' % (tbin,offset)
                    tbin += 1

            hist[binnum] += 1
            binnum_old = binnum
        #print 'rev:',revind[binnum]

        i += 1
        offset += 1

    if dorev:
        # Fill in the last ones
        tbin = binnum_old + 1
        while tbin <= nbin:
            revind[tbin] = revind.size
            tbin += 1





def histogram2d(x, y, 
                nx=None, 
                ny=None, 
                xbin=None, 
                ybin=None, 
                xmin=None, 
                xmax=None, 
                ymin=None, 
                ymax=None, 
                rev=False):

    if x.size != y.size:
        raise ValueError("x and y must be the same size")

    # binsizes will take precedence
    dobinsizes=False
    donbin=False
    if (xbin != None) and (ybin != None):
        dobinsizes=True
    elif (nx != None) or (ny != None):
        donbin=True
    else:
        raise ValueError("Enter either nx,ny or xbin,ybin")

    if xmin is None:
        xmin=x.min()
    if xmax is None:
        xmax=x.max()
    if ymin is None:
        ymin=y.min()
    if ymax is None:
        ymax=y.max()

    w, = numpy.where( (x >= xmin) & (x <= xmax) & (y >= ymin) & (y <= ymax) )
    if w.size == 0:
        raise ValueError("No data in specified min/max range\n")

    if dobinsizes:
        # determine nx,ny from binsizes
        nx = numpy.int64(  (xmax-xmin)/xbin ) + 1
        ny = numpy.int64(  (ymax-ymin)/ybin ) + 1

    xind=numpy.floor((x[w]-xmin)*(nx/(xmax-xmin)))
    yind=numpy.floor((y[w]-ymin)*(ny/(ymax-ymin)))

    ind=xind+nx*yind


    hist=histogram(ind, min=0l, max=nx*ny-1, rev=rev)
    if rev:
        hist, revind=hist
        hist = hist.reshape(nx,ny)
        return hist, revind
    else:
        hist = hist.reshape(nx,ny)
        return hist


def wmom(arrin, weightsin, inputmean=None, calcerr=False, sdev=False):
    """
    NAME:
      wmom()
      
    PURPOSE:
      Calculate the weighted mean, error, and optionally standard deviation
      of an input array.

    CALLING SEQUENCE:
     wmean,werr = wmom(arr, weights, inputmean=None, calcerr=False, sdev=False)
    
    INPUTS:
      arr: A numpy array or a sequence that can be converted.
      weights: A set of weights for each elements in array.
    OPTIONAL INPUTS:
      inputmean: An input mean value, around which them mean is calculated.
      calcerr=False: Calculate the weighted error.  By default the error
        is calculated as 1/sqrt( weights.sum() ).  In this case it is
        calculated as sqrt( (w**2 * (arr-mean)**2).sum() )/weights.sum()
      sdev=False: If True, also return the weighted standard deviation 
        as a third element in the tuple.

    OUTPUTS:
      wmean, werr: A tuple of the weighted mean and error. If sdev=True the
         tuple will also contain sdev: wmean,werr,wsdev

    REVISION HISTORY:
      Converted from IDL: 2006-10-23. Erin Sheldon, NYU

   """
    from numpy import float64
    
    # no copy made if they are already arrays
    arr = numpy.array(arrin, ndmin=1, copy=False)
    weights = numpy.array(weightsin, ndmin=1, copy=False)
    
    # Weights is forced to be type double. All resulting calculations
    # will also be double
    if weights.dtype != float64:
        weights = numpy.array(weights, dtype=float64)
  
    wtot = weights.sum()
        
    # user has input a mean value
    if inputmean is None:
        wmean = ( weights*arr ).sum()/wtot
    else:
        wmean=float(inputmean)

    # how should error be calculated?
    if calcerr:
        werr2 = ( weights**2 * (arr-wmean)**2 ).sum()
        werr = numpy.sqrt( werr2 )/wtot
    else:
        werr = 1.0/numpy.sqrt(wtot)

    # should output include the weighted standard deviation?
    if sdev:
        wvar = ( weights*(arr-wmean)**2 ).sum()/wtot
        wsdev = numpy.sqrt(wvar)
        return wmean,werr,wsdev
    else:
        return wmean,werr




def sigma_clip(arrin, niter=4, nsig=4, extra={}, verbose=False):
    """
    NAME:
      sigma_clip()
      
    PURPOSE:
      Calculate the mean/stdev of an array with sigma clipping. Iterate
      niter times, removing elements that are outside nsig, and recalculating
      mean/stdev.

    CALLING SEQUENCE:
      mean,stdev = sigma_clip(arr, niter=4, nsig=4, extra={})
    
    INPUTS:
      arr: A numpy array or a sequence that can be converted.

    OPTIONAL INPUTS:
      niter: number of iterations, defaults to 4
      nsig: number of sigma, defaults to 4

    OUTPUTS:
      mean,stdev: A tuple containing mean and standard deviation.
    OPTIONAL OUTPUTS
      extra={}: Dictionary containing the array of used indices in
         extra['index']

    REVISION HISTORY:
      Converted from IDL: 2006-10-23. Erin Sheldon, NYU

   """
    arr = numpy.array(arrin, ndmin=1, copy=False)

    index = numpy.arange( arr.size )

    for i in numpy.arange(niter):
        m = arr[index].mean()
        s = arr[index].std()

        if verbose:
            sys.stdout.write('iter %s\tnuse: %s\tmean %s\tstdev %s\n' % \
                (i+1, index.size,m,s))

        clip = nsig*s

        w, = numpy.where( (numpy.abs(arr[index]) - m) < clip )

        if w.size == 0:
            sys.stderr.write("nsig too small. Everything clipped on iteration %d" % i+1)
            return m,s


        index = index[w]

    # Calculate final stats
    amean = arr[index].mean()
    asig = arr[index].std()

    extra['index'] = index
    return amean, asig
     

