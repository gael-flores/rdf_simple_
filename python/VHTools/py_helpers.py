import ROOT
import numpy as np
import glob
from contextlib import contextmanager
from common.plotter import *
import ast

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~plotting with matplotlib
def color(kColor):
    """
    Retrieve the normalized RGB components of a ROOT color.

    Parameters:
    ----------
    kColor : int
        ROOT color index (e.g., ROOT.kRed, ROOT.kBlue+2, etc.).

    Returns:
    -------
    tuple of float
        A tuple (r, g, b) where each value is in the range [0.0, 1.0].

    Raises:
    ------
    ValueError
        If the color index does not correspond to a valid ROOT TColor.

    Example:
    --------
    >>> color(ROOT.kGreen+2)
    (0.0, 0.6, 0.0)
    """
    color = ROOT.gROOT.GetColor(kColor)
    if not color:
        raise ValueError(f"Invalid ROOT color index: {kColor}")

    r = color.GetRed()
    g = color.GetGreen()
    b = color.GetBlue()

    return (r, g, b)

def getPoisson_arrays(h, scale=1):
    """
    Extract Poisson statistics (Garwood interval) from a ROOT TH1 histogram 
    and return arrays for plotting with asymmetric error bars.

    Parameters:
    ----------
    h : ROOT.TH1
        Input histogram with Poisson-distributed bin contents.

    scale : float, optional
        Scaling factor to apply to bin contents and errors (default is 1).
        Useful when normalizing or converting to different units (e.g., cross-section).

    Returns:
    -------
    tuple of np.ndarray
        A 5-element tuple:
        - bin_centers: array of bin centers (x positions)
        - bin_contents: array of (possibly scaled) bin contents (y values)
        - lower_errors: array of (possibly scaled) lower Poisson errors
        - upper_errors: array of (possibly scaled) upper Poisson errors
        - bin_widths: array of bin widths (useful for bar-style plotting)

    Notes:
    -----
    - Poisson errors are computed using the Garwood method via `fetchError(q, N)`.
    - Bins with zero content will still return asymmetric error bars.
    - Designed for use with plotting libraries like Matplotlib or mplhep.

    Example:
    --------
    >>> x, y, yerr_low, yerr_up, width = getPoisson_arrays(hist)
    >>> plt.errorbar(x, y, yerr=[yerr_low, yerr_up], fmt='o')
    """
    q = (1 - 0.6827) / 2.0  # 68.27% confidence interval (1σ)

    bin_centers  = []
    bin_contents = []
    bin_widths   = []
    lower_errors = []
    upper_errors = []

    for i in range(1, h.GetNbinsX() + 1):
        thresh = h.GetXaxis().GetBinCenter(i)
        N = h.GetBinContent(i)
        bin_width = h.GetXaxis().GetBinWidth(i)

        error_low, error_up = fetchError(q, N)

        bin_centers.append(thresh)
        bin_contents.append(N)
        bin_widths.append(bin_width)
        lower_errors.append(N - error_low)
        upper_errors.append(error_up - N)

    return (
        np.array(bin_centers),
        scale * np.array(bin_contents),
        scale * np.array(lower_errors),
        scale * np.array(upper_errors),
        np.array(bin_widths)
    )

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ROOT histograms and plotting
def canvas_partition2(canvas, lmargin=0.1, rmargin=0., tmargin=0.07, bmargin=0.10, frame_line_width=3):

    """
    Partition a ROOT canvas into 1 row and 4 vertical pads (columns),
    with a wider leftmost pad to accommodate a Y-axis label.

    Parameters:
    ----------
    canvas : ROOT.TCanvas
        The canvas to partition. It should be empty when passed in.

    lmargin : float, optional
        Left margin for the leftmost pad (default is 0.1).

    rmargin : float, optional
        Right margin for the rightmost pad (default is 0.0).

    tmargin : float, optional
        Top margin for all pads (default is 0.07).

    bmargin : float, optional
        Bottom margin for all pads (default is 0.10).

    frame_line_width : int, optional
        Width of the pad frame lines (default is 3).

    Returns:
    -------
    dict
        Dictionary of pads indexed by (i, 0) where `i` is the column index (0 to 3).

    Notes:
    -----
    - Leftmost pad is wider to allow for axis labels.
    - Internal horizontal spacing is handled by setting individual widths.
    - You can access a specific pad via: `pads[(i, 0)]`
    - The pads are numbered from 1 to 4 internally via `SetNumber`.

    Example:
    --------
    >>> c = ROOT.TCanvas("c", "", 1200, 400)
    >>> pads = canvas_partition(c)
    >>> pads[(0, 0)].cd()
    >>> hist.Draw()
    >>> c.Update()
    """
    pads = {}
    total_pads = 4
    # Define relative width of each pad
    left_pad_frac = 0.25
    mid_left_frac = 0.24
    mid_right_frac = 0.24
    right_pad_frac = 0.27
    other_pad_frac = 0.24

    x_start = 0.0
    for i in range(total_pads):
        # Set width depending on position
        if i == 3:
            width = right_pad_frac
        elif i == 0:
            width = left_pad_frac
        else:
            width = other_pad_frac

        xlow = x_start
        xup = x_start + width
        x_start = xup

        pad = ROOT.TPad(f"pad_{i}_0", "", xlow, 0.0, xup, 1.0)
        pad.SetNumber(i + 1)

        # Set margins: left pad gets left margin, right pad gets right margin
        pad.SetLeftMargin(lmargin if i == 0 else 0)
        pad.SetRightMargin(rmargin if i == total_pads - 1 else 0)
        pad.SetTopMargin(tmargin)
        pad.SetBottomMargin(bmargin)

        # Frame styling
        pad.SetFrameLineWidth(frame_line_width)
        pad.SetFrameLineColor(ROOT.kBlack)

        pad.Draw()
        pads[(i, 0)] = pad

    return pads

def x_axis(col,binsM=1):
    """
    Create and customize a TGaxis for labeling a custom x-axis in prefit and postfit grid plots.

    This function is designed to label only even-numbered positions in the 11*binsM subdivision.

    Parameters:
    ----------
    binsM : int
        Number of mass bins or base subdivisions. The total number of label slots 
        will be 11 * binsM. Only a subset of those will be labeled to reduce clutter.

    col : int
        Selector for axis labeling style:
        - If col == 3, the axis will have the title "L_{xy} [cm]".
        - Otherwise, no title is applied.

    Returns:
    -------
    ROOT.TGaxis
        The constructed and drawn TGaxis object.

    Notes:
    -----
    - Labels are placed at every 20th bin (i.e., every 2nd of 11 minor divisions).
    - Axis title and styling are hardcoded for visual appearance.
    - The axis is drawn along the bottom edge of the current pad (gPad).
    - It assumes the user has already set up a canvas and drawn something else,
      so this axis overlays an existing frame.

    Example:
    --------
    >>> c = ROOT.TCanvas()
    >>> h.Draw()
    >>> x_axis(10, 3)
    """
    # Create axis from lower left to lower right of pad
    ax = ROOT.TGaxis(
        ROOT.gPad.GetUxmin(), ROOT.gPad.GetUymin(),
        ROOT.gPad.GetUxmax(), ROOT.gPad.GetUymin(),
        0, binsM * 11, binsM * 11, "+"
    )

    # Customize labels: only label every other major division (i.e., even multiples)
    for n in range(0, 11 * binsM):
        if binsM == 1: 
            if (n%11)%1 == 0:
                ax.ChangeLabel(n+1, -1, 0.025, -1, -1, -1, "{}".format((n*10)%110-10))
            else:
                ax.ChangeLabel(n+1, -1, 0.0, -1, -1, -1, "")    
            if col == 1:
                ax.SetTitle("L_{xy} [cm]")
                ax.ChangeLabel(11+1, -1, 0.025, -1, -1, -1, "100")
            else:
                ax.ChangeLabel(11+1, -1, 0.025, -1, -1, -1, "-10")
        else:    
            if (n % 11) % 2 == 0:
                ax.ChangeLabel(n + 1, -1, 0.023, -1, -1, -1, "{}".format((n * 10) % 110 - 10))
            else:
                ax.ChangeLabel(n + 1, -1, 0.0, -1, -1, -1, "")  # hide label
            if col == 1:
                ax.SetTitle("L_{xy} [cm]")
                ax.ChangeLabel(11*n +1, -1, 0.0, -1, -1, -1, "")
                ax.ChangeLabel(1, -1, 0.023, -1, -1, -1, "-10")
            else:
                ax.ChangeLabel(11*n + 1, -1, 0.023, -1, -1, -1, "-10")

    # Styling
    ax.SetLabelSize(0.1)
    ax.SetLabelOffset(-0.055)
    #ax.SetLabelFont(2)

    # Draw the axis on top of the current pad
    ax.Draw("same")
    return ax

def m_range(binsM,m):
    """
    Draw a custom TGaxis along the top of the current ROOT pad, labeling 
    mass bins in the form: "4.00 ≤ m(γ,γ) < m + 5 GeV".

    Parameters:
    ----------
    binsM : int
        Number of mass bins to divide the range into.

    m : float
        Total mass range width. The range is assumed to start at 4.0 GeV 
        and end at 5.0 + m GeV.

    Returns:
    -------
    ROOT.TGaxis
        The TGaxis object drawn at the top of the current canvas.

    Notes:
    -----
    - The mass range is from 4.0 to (5.0 + m) GeV.
    - The axis overlays the top of the current pad (gPad).
    - Labels are centered in each mass bin and use compact formatting.
    - Label size and offset are set for visual clarity.

    Example:
    --------
    >>> m_range(10, 5.0)  # 10 bins from 4.0 to 10.0 GeV
    """
    ax_m = ROOT.TGaxis(
        ROOT.gPad.GetUxmin(), ROOT.gPad.GetUymax(),
        ROOT.gPad.GetUxmax(), ROOT.gPad.GetUymax(),
        0, binsM, binsM, "-")
    
    for n in range(binsM):
        ax_m.ChangeLabel(n+1, -1, 0.02, -1, -1, -1, f"{4+(m+1.)/binsM*n:.2f} #leq m(#gamma,#gamma) < {4+(m+1.)/binsM*(n+1):.2f} GeV")
    ax_m.CenterLabels()
    ax_m.SetLabelOffset(-0.08)
    ax_m.Draw("same")
    return ax_m

def get_max_y_with_error(graph):
    """
    Compute the maximum Y value from a ROOT.TGraphAsymmErrors object.

    Parameters:
    ----------
    graph : ROOT.TGraphAsymmErrors
        Input graph with asymmetric error bars.

    Returns:
    -------
    float
        Maximum value of (Y + upper error) across all points in the graph.

    Example:
    --------
    >>> max_y = get_max_y_with_error(g)
    >>> h.SetMaximum(1.1 * max_y)  # add padding when drawing
    """
    try:
        n = graph.GetN()
        max_y = float('-inf')
        for i in range(n):
            y = graph.GetPointY(i)
            y_err_up = graph.GetErrorYhigh(i)
            y_total = y + y_err_up
            if y_total > max_y:
                max_y = y_total
    except AttributeError:
        max_y = graph.GetMaximum() + graph.GetBinErrorUp(graph.GetMaximumBin())
    return max_y

def getPoisson(h,scale=1):

    """
    Convert a ROOT TH1 histogram into a TGraphAsymmErrors object with 
    asymmetric Poisson confidence intervals (Garwood intervals) as error bars.
    Optional scaling as second argument


    Parameters:
    ----------
    h : ROOT.TH1
        Input histogram. Bin contents are assumed to follow Poisson statistics.
    scale: float
        Value to scale histogram by.

    Returns:
    -------
    ROOT.TGraphAsymmErrors
        A TGraphAsymmErrors object with:
        - X values = bin centers
        - Y values = bin contents
        - Asymmetric Y errors = Poisson (Garwood) confidence intervals
        - No X error bars

    Notes:
    -----
    - Uses a 68.27% confidence level by default (1σ).
    - Errors are computed even for bins with zero entries (lower error = 0, upper error > 0).
    - Points with zero content are included in the graph (not skipped).
    - `fetchError()` is assumed to return [lower, upper] bounds of the Poisson interval.

    Example:
    --------
    >>> h = ROOT.TH1F("h", "", 10, 0, 10)
    >>> h.Fill(1)
    >>> g = getPoisson(h)
    >>> g.Draw("AP")
    """
    q = (1-0.6827)/2.0
    gRate = ROOT.TGraphAsymmErrors()
    n = 0
    for i in range(1, h.GetNbinsX()+1):
        thresh = h.GetBinCenter(i)
        N = h.GetBinContent(i)
        #if N == 0:
        #    continue
        gRate.SetPoint(n, thresh, N * scale)
        error = fetchError(q, N)
        print(error)
        if scale == 1:
            gRate.SetPointError(n, h.GetBinWidth(i)/2.0, h.GetBinWidth(i)/2.0, (N-error[0]), (error[1]-N))
        else:
            gRate.SetPointError(n, h.GetBinWidth(i)/2.0, h.GetBinWidth(i)/2.0, scale*(N-error[0]), scale*(error[1]-N))
        n+=1
    return gRate

def getPoisson_from_DF(h,df,label_low,label_high):
    gRate = ROOT.TGraphAsymmErrors()
    n = 0
    for i in range(1, h.GetNbinsX()+1):
        err_low = df.loc[n][label_low]
        err_up  = df.loc[n][label_high]
        thresh = h.GetBinCenter(i)
        N = h.GetBinContent(i)
        gRate.SetPoint(n, thresh, N)
        gRate.SetPointError(n, h.GetBinWidth(i)/2.0, h.GetBinWidth(i)/2.0, (err_low), (err_up))
        n+=1
    return gRate

def getPoisson_fromtxt(h,v,l,m,era,sfx):
    gRate = ROOT.TGraphAsymmErrors()
    n = 0
    for i in range(1, h.GetNbinsX()+1):
        thresh = h.GetBinCenter(i)
        N = h.GetBinContent(i)
        gRate.SetPoint(n, thresh, N) #weighted value
        error = fetchError_fromtxt(v,l,m,era,n,sfx)
        gRate.SetPointError(n,h.GetBinWidth(i)/2.0, h.GetBinWidth(i)/2.0, (error[0]), (error[1]))
        
        n+=1
    return gRate

def fillZeroBins(hist, thresh):
    """
    Replace zero or negative bin contents in a ROOT histogram with a threshold value.

    This is used when:
    - Avoiding issues with logarithmic plotting (where log(0) is undefined).
    - Ensuring non-zero content for statistical operations or plotting.

    Parameters:
    ----------
    hist : ROOT.TH1
        The input histogram whose bins will be modified in-place.

    thresh : float
        The value to assign to bins with content <= 0.

    Returns:
    -------
    ROOT.TH1
        The modified histogram (same object, modified in-place).

    Example:
    --------
    >>> h = fillZeroBins(h, 1e-3)
    """
    for i in range(1, hist.GetNbinsX() + 1):
        if hist.GetBinContent(i) <= 0:
            hist.SetBinContent(i, thresh)
    return hist

def canvas_partition(canvas, lmargin=0.1, rmargin=0., tmargin=0.07, bmargin=0.10, frame_line_width=2):
    """Partitions canvas into 1 row and 4 columns, with wider leftmost pad for Y-axis."""
    pads = {}
    
    total_pads = 2
    
    left_pad_frac = 0.48
    right_pad_frac = 0.48
    other_pad_frac = 0.5
    #other_pad_frac = (1.0 - left_pad_frac) / (total_pads - 1)

    x_start = 0.0
    for i in range(total_pads):
        if i == 1:
            width = right_pad_frac
        elif i == 0:
            width = left_pad_frac
        else:
            width = other_pad_frac

        xlow = x_start
        xup = x_start + width
        x_start = xup

        pad = ROOT.TPad(f"pad_{i}_0", "", xlow, 0.0, xup, 1.0)
        pad.SetNumber(i + 1)

        # Margin logic
        pad.SetLeftMargin(0.1 if i == 0 else 0)
        #pad.SetRightMargin(0.01 if i == total_pads - 1 else 0)
        #for right axis label:
        pad.SetRightMargin(0.13 if i == total_pads - 1 else 0)
        pad.SetTopMargin(tmargin)
        pad.SetBottomMargin(bmargin)

        pad.SetFrameLineWidth(frame_line_width)
        pad.SetFrameLineColor(ROOT.kBlack)

        pad.Draw()
        pads[(i, 0)] = pad

    return pads

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~calculations
def fetchError(q, n):
    """
    Compute the lower and upper bounds of a confidence interval using 
    the inverse chi-squared distribution (Garwood interval).

    This function is commonly used to determine Poisson confidence intervals 
    for an observed number of events `n` at a given confidence level.

    Parameters:
    ----------
    q : float
        Quantile value for the desired confidence interval (e.g., 0.16 for 68% CL).
        Note: For 68% central intervals, use q = alpha/2 = 0.16.

    n : int
        Number of observed events (must be >= 0).

    Returns:
    -------
    list of float
        A list [lower, upper] representing the lower and upper bounds 
        of the confidence interval.

    Notes:
    -----
    - If n = 0, the lower bound is set to 0.
    - Uses the central quantile of the chi-squared distribution for interval construction.
    - This is equivalent to the Garwood construction of a Poisson confidence interval.

    Example:
    --------
    >>> fetchError(0.16, 3)
    [0.584, 7.255]  # (values are approximate)
    """
    l = 0
    if n!=0:
        l = ROOT.Math.chisquared_quantile_c(1.0-q, 2.0*n)/2.0
    u = ROOT.Math.chisquared_quantile_c(q, 2.0*n+2)/2.0
    return [l,u]

def fetchError_fromtxt(v,l,m,era,n,sfx):
    #with open(f'/uscms/home/gfavila/nobackup/rdf_simple_/{v}H_{l}_m{m}_{era}_WGammaplusTTJets.txt') as f:
    with open(f'/uscms/home/gfavila/nobackup/rdf_simple_/{v}H_{l}_m{m}_{era}_{sfx}.txt') as f:
        data = f.read()
    intervals = ast.literal_eval(data)
    return intervals[v][l][m][n]

def chi_squared_graphs(g1, g2):
    """
    Compute the chi-squared between two TGraphAsymmErrors.

    Parameters:
        g1, g2 (ROOT.TGraphAsymmErrors): The graphs to compare.
            Typically g1 is data and g2 is model.

    Returns:
        chi2 (float): Total chi²
        ndf (int): Number of points used
    """
    n1 = g1.GetN()
    n2 = g2.GetN()

    if n1 != n2:
        raise ValueError("Graphs must have the same number of points")

    chi2 = 0.0
    ndf = 0

    for i in range(n1):
        x1, y1 = np.array([0.0]), np.array([0.0])
        x2, y2 = np.array([0.0]), np.array([0.0])
        g1.GetPoint(i, x1, y1)
        g2.GetPoint(i, x2, y2)

        # You could add a check here to compare x1 and x2, if needed
        if x1 != x2:
            print(f"Warning: x mismatch at index {i}: x1 = {x1}, x2 = {x2}")

        # Combine asymmetric errors in quadrature
        e1_low = g1.GetErrorYlow(i)
        e1_high = g1.GetErrorYhigh(i)
        

        # Choose upper or lower error depending on which side the difference is
        err1 = e1_high #if y1 > y2 else e1_low
        

        total_err = err1
        if total_err == 0:
            continue  # skip to avoid divide by zero

        chi2 += ((y1 - y2)/total_err)**2 
        ndf += 1

    reduced_chi2 = chi2 / ndf if ndf > 0 else float("inf")
    
    return chi2[0], ndf

def chi_squared(hist1, hist2, ignore_empty_bins=True):
    """
    Calculate the chi-squared value between two ROOT.TH1D histograms.

    Parameters:
        hist1 (ROOT.TH1D): First histogram (observed).
        hist2 (ROOT.TH1D): Second histogram (expected - background).
        ignore_empty_bins (bool): If True, skip bins where both histograms have 0 content.

    Returns:
        chi2 (float): The chi-squared value.
        ndf (int): Number of degrees of freedom (number of bins used).
    """
    if hist1.GetNbinsX() != hist2.GetNbinsX():
        raise ValueError("Histograms must have the same number of bins")

    chi2 = 0.0
    
    ndf = 0

    for i in range(1, hist1.GetNbinsX() + 1):  # bin numbering starts at 1 in ROOT
        obs = hist1.GetBinContent(i)
        exp = hist2.GetBinContent(i)
        e1 =  hist1.GetBinError(i)
        if ignore_empty_bins and obs == 0 and exp == 0:
            continue
        
        if e1 == 0:
            continue  # skip bin if total error is zero to avoid division by zero

        chi2 += ((obs - exp)/e1)**2
        ndf += 1

    return chi2, ndf

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~file management
@contextmanager
def open_root_file(filename, mode="READ"):
    """
    Context manager for safely opening and closing a ROOT TFile.

    This ensures that the ROOT file is automatically closed when the
    context exits, even in the event of an exception.

    Parameters:
    ----------
    filename : str
        Path to the ROOT file to open.

    mode : str, optional
        Mode in which to open the file (default is "READ").
        Typical options include:
        - "READ"   : Open existing file for reading
        - "UPDATE" : Open existing file for reading/writing
        - "RECREATE": Create new file, overwriting if exists
        - "NEW"    : Create new file, fail if exists

    Yields:
    ------
    ROOT.TFile
        The opened ROOT file object.

    Example:
    --------
    >>> with open_root_file("data.root") as f:
    ...     h = f.Get("myHistogram")
    ...     h.Draw()
    """
    file = ROOT.TFile.Open(filename, mode)
    try:
        yield file
    finally:
        if file and file.IsOpen():
            file.Close()

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~file management
