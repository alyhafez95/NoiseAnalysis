# Noise Analysis
* this program was designed to be flexable and be used as a way to analyze noise from all types of analytical assays.
* doubleclick on NoiseAnalysis.exe file to start.
* a windows file explorer will pop up > select desired .csv or a .arw file to with your sequence blanks (or whatever you want to analyze).
* the columns' format must be as such: x(1), y(1), x(2) , y(2), ... x(n), y(n)...etc.
* you may select a .csv or a .arw file with only one blank at a time if you prefer: x(1), y(1).
* you may also selct multiple files at a time.
* answer the prompt questions to guide the program.
    - the following are acceptable inputs for a True/Yes answer (case insenstive):
            ['1', 't', 'tr', 'tru', 'true', 'truee', 'y', 'ye', 'yes', 'yess', 'yeah', 'yu', 'yup', 'yupp', 'sure', 'certainly', 'yay']
    - any other input will be regarded as a False/No (including no input).
* the script will take on the input data by sets of 2 columns at a time x(n), y(n) and do the following steps:
    <!-- 1. (removed) if the number of data points is lower the 25,000 the script will over sample the data via interpolation till that is no longer the case.
        * this is to establish some consistancy in the number of data points -->
    1. (optional) the script will then run its "automatic peak detection":
        * this is done via a smoothing algorithom that runs along the xy data: 
            SMA_i[i] = sum(y[i-range:i+range])/len(y[i-range:i+range])
        * where the smoothed curve (SMA_i[i]) and the actual data (y[i]) are decoupled by a number greater than a pre-selected threshold (a multiple of the data's standard deviation) y[i] is then characterized as a signal.
        * any false negatives points are corrected by conforming them to their nearest n neighboors
        * keep in mind this will also characterize negative peaks as signals.
        * points characerized as signals will be removed, what remains should be just the baseline.
    2. (optional) the script will preform a polynomial fit of the baseline (default deg=4)
        * (optional) data points further away from the fit will be eliminated.
        * data points will then be flattened by subtracting the polynomial fit from the left over baseline
    3. (optional) a fourier transformation will be done on the remaining data, creating a full fourier series representation of the data. This is a way to eliminate some if not all the pattern oscillations in a baseline (use sparingly).
        * the baseling then will be approximated using a user selected number of harmonics; and the baseline will once again be flattened via subtraction of the reconstructed approximation of the baseline shape
    4. In the final step, the script will evaluate the left over baseline (y2) and make the following calculations:
        noise = y2.std()
        LoD = 3*y2.std()
        LoQ = 10*y2.std()

        Note: The script will always center the data around zero; and when the mean = 0, RMS = Standard Deviation.
              LoD and LoQ RMS multiples are per USP guidlines.

    # Output:
    - _rms_noise.csv file will appear in the same directory as the data you selected, as well as .png files with reports on each individual set of x,y data.
    - if the reports indicate a bad replicate you should remove that replicate from your final LoD/LoQ average.
    - you may also elect to or rerun the script with adjustments to the parameters.
    - it is recommended however to find a set of parameters that work well for your assay and stick to those parameters as much as possible to maintain longterm consistancy.
