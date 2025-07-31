# Function for combined normalised plot and normalised ratio plot 
import matplotlib as mpl 
import matplotlib.pyplot as plt


def errors_on_ratio_of_normalised(list_of_hists,index_of_denom=0):
    
    # Normalise the histograms to unity
    norm_hist = lambda x: x/x.sum()
    list_of_norm_hists = list(map(norm_hist, list_of_hists))
    
    # Calculate the errors on the normalised histograms
    norm_hist_error = lambda unnorm_hist:  unnorm_hist.variances()**0.5/unnorm_hist.sum()
    errors = list(map(norm_hist_error, list_of_hists))
    
    # Calculate the ratio of the normalised histograms
    # and the errors on the ratio
    # The ratio is calculated as ratio_hist = hist_i / hist_denom
    # The errors on the ratio are calculated as:
    # ratio_hist_errors = ratio_hist * ((errors_denom/hist_denom)**2 + (errors_i/hist_i)**2)**0.5
    ratio_hist_list =[]
    ratio_error_list = []
    for i in range(len(list_of_norm_hists)):
        ratio_hist = list_of_norm_hists[i]/list_of_norm_hists[index_of_denom]
        
        numerator_term = errors[i]/list_of_norm_hists[i].view()
        denominator_term = errors[index_of_denom]/list_of_norm_hists[index_of_denom].view()
        
        ratio_hist_errors = ratio_hist.view() * (numerator_term**2 + denominator_term**2 )**0.5
    
        ratio_hist_list.append(ratio_hist)
        ratio_error_list.append(ratio_hist_errors)
    return ratio_hist_list, ratio_error_list
    


def normalised_ratio_plot(list_of_hists):

    hists,errors = errors_on_ratio_of_normalised(list_of_hists, index_of_denom=0)
    colors = plt.rcParams['axes.prop_cycle'].by_key()['color']

    fig, ax = plt.subplots()
    hep.histplot(hists, label=["SM", "cQQ8 up", "cQQ8 down"])
    for hist,err,c in zip(hists, errors,colors):
        plt.errorbar(hist.axes[0].centers,
                    hist.values(),
                    err, fmt='none',color=c)
        
    return fig,ax



def combined_normalised_plot_and_ratio(list_of_hists):
    
    """
    Plotting script for normalised histograms and normalised ratio plot.
    """

    mpl.rcParams.update(mpl.rcParamsDefault)
    plt.style.use('default')
    
    # turn on tex style 

    plt.rcParams.update({
        "font.family": "serif",
        "font.serif": ["Palatino", "Times New Roman", "Times", "CMU Serif"],
        "mathtext.fontset": "cm",  # Use Computer Modern fonts for math text (built-in)
        "mathtext.rm": "serif",    # Use serif fonts for regular math text
    })


    fig, ax = plt.subplots(
    2, 1,           # 2 rows, 1 column
    figsize=(6, 6), # Width x Height in inches
    gridspec_kw={'height_ratios': [2.5, 1]}, # Top:Bottom height ratio
    dpi=200,
    sharex=True     # Share x-axis
)
    plt.subplots_adjust(hspace=0.05)
    colors = plt.rcParams['axes.prop_cycle'].by_key()['color']

    # Top plot - histograms normalised to unity
    norm_hist = lambda x: x/x.sum()
    list_of_norm_hists = list(map(norm_hist, list_of_hists))
    
    for i in range(len(list_of_norm_hists)):
        ax[0].errorbar(list_of_norm_hists[i].axes[0].centers,
                     list_of_norm_hists[i].values(),
                     list_of_hists[i].variances()**0.5/list_of_hists[i].sum(), fmt='none', color=colors[i])
    
    hep.histplot(list_of_norm_hists,ax=ax[0])
    ax[0].legend(["SM", "cQQ8 up", "cQQ8 down"])
    ax[0].set_xlabel("")
    ax[0].set_ylabel("Normalised events")
    # Remove x-axis label from top plot 
    # Normalised ratio plot
    hists, errors = errors_on_ratio_of_normalised(list_of_hists, index_of_denom=0)
    hep.histplot(hists, label=["SM", "cQQ8 up", "cQQ8 down"], ax=ax[1])
    for hist, err, c in zip(hists, errors, colors):
        ax[1].errorbar(hist.axes[0].centers,
                     hist.values(),
                     err, fmt='none',color=c)
    
    ax[1].set_xlabel("Random top pair combined mass [GeV]")
    ax[1].set_ylabel("Ratio to SM")
    
    ax[0].tick_params(direction='in')
    ax[1].tick_params(direction='in')
    
    return fig, ax