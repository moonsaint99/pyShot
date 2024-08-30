# This script lets you use an interactive plot where you can use your cursor to pick peaks
# in the seismic signal, and the script will snap to the nearest peak/trough.
# You don't have to pick every single trace, just the ones that look like they
# have the arrival you need to pick.
# It will save your picks to a csv file.
# The csv file will include pick time, pick amplitude, standard deviation in a window around your pick,
# trace offset, and trace number.
#

#from ..load.load_segy import loadshots

import matplotlib as mpl
# mpl.use('macosx')
import matplotlib.pyplot as plt
import obspy as op
import numpy as np
import scipy as sp
import os

# # Load all .su files in the directory of choice into a dictionary of ObsPy Stream objects.
# lithic_smallstack_streams = loadshots("lithic_smallstack_streams/abs/")

# Helper functions
def add_pick(pick, pickdict):
    pickdict[pick[0]] = (pick[1], pick[2], pick[3])
    pickdict = dict(sorted(pickdict.items()))
    print(pickdict)

def clip_stream(stream_original, clipval):
    stream = stream_original.copy()
    for trace in stream:
        # trace.data = np.clip(trace.data, -1*clipval, clipval)
        trace.data = np.where((trace.data >= -1*clipval) & (trace.data <= clipval), trace.data, 0)
    return stream

# Define a function to plot a Stream object and let you pick peaks in the traces.
# n and m to increase or decrease the gain.
# You can click on a trace to pick the nearest peak or trough.

class PickPlotParams:
    """
    A class for storing plotting parameters for the PickPlotter class.
    """

class PickPlotter:
    """
    A class for interactively picking seismic arrivals from a loaded stream object.

    """

    def __init__(self, stream_input, filename="Record Section", outfile="picks.csv"):
        """
        Initializes the PickPlotter class with stream data, filename, and output filename.

        Args:
            stream_input: The target obspy stream object.
            filename: The name of the seismic record.
            outfile: The filename for the output pick csv file.
        """
        # Stream object, both visual and original
        self.stream_visual = stream_input.copy()
        self.stream_original = stream_input.copy()
        # Name and directory of output pick file
        self.filename = filename
        self.outfile_name = os.path.basename(outfile)
        self.outfile_dir = os.path.dirname(outfile)

        # Plotting parameters
        self.data_delta = self.stream_original[0].stats.delta
        self.scale = 1

        # Clip value for traces
        trace_maxvals = []
        for trace in self.stream_visual:
            trace_maxvals.append(np.max(np.abs(trace.data)))
        maxamp = np.max(trace_maxvals)
        self.clip = maxamp

        # Peak information
        self.peaks = []
        self.peaktimes = []
        for trace in self.stream_original:
            peakindices = sp.signal.find_peaks(trace.data)[0]
            troughindices = sp.signal.find_peaks(-trace.data)[0]
            supremaindices = np.concatenate((peakindices, troughindices))
            self.peaks.append(supremaindices)
            suprematimes_array = trace.stats.delta * supremaindices
            self.peaktimes.append(suprematimes_array.tolist())

        self.pickdictPlot1 = None
        self.pickdictPlot2 = None

        #initial plotting --- all AI stuff below this

        self.ax = plt.gca()
        self.ax.set_title(self.filename)
        self.ax.set_xlabel('Offset (km)')
        self.ax.set_ylabel('Time (s)')
        self.ax.set_ylim(np.require(self.stream_visual[0].stats.endtime, dtype='float32'), 0)

    def load_data(self):
        """
        Precomputes peak information for all traces in the stream.
        """
        self.data_delta = self.stream[0].stats.delta
        peaks = []
        peaktimes = []
        for trace in self.stream_original:
            peakindices = sp.find_peaks(trace.data)[0]
            troughindices = sp.find_peaks(-trace.data)[0]
            supremaindices = np.concatenate((peakindices, troughindices))
            peaks.append(supremaindices)
            suprematimes_array = trace.stats.delta * supremaindices
            peaktimes.append(suprematimes_array.tolist())
        self.peaks = peaks
        self.peaktimes = peaktimes

    def plot_data(self):
        """
        Plots the seismic section data and initializes interactive elements.
        """
        self.fig, self.ax = plt.subplots(figsize=(12, 8))
        self.stream.plot(
            type="section",
            scale=self.scale,
            fig=self.fig,
            time_down=True,
            fillcolors=("blue", "red"),
            color="black",
            size=(600, 800),
        )
        self.ax.set_title(self.filename)
        self.ax.set_xlabel("Offset (km)")
        self.ax.set_ylabel("Time (s)")
        self.ax.set_ylim(np.require(self.stream[0].stats.endtime, dtype="float32"), 0)

        offsetlist = [trace.stats.distance / 1000 for trace in self.stream_original]
        self.point, = self.ax.plot([], [], "_", markersize=25, color="green")
        self.pickInfo = self.ax.text(
            0.02,
            0.02,
            "",
            transform=self.ax.transAxes,
            verticalalignment="top",
            bbox=dict(boxstyle="round", facecolor="wheat", alpha=1),
        )

        self.fig.canvas.mpl_connect("key_press_event", self.handle_key_press)
        self.fig.canvas.mpl_connect("button_press_event", self.handle_mouse_click)
        self.fig.canvas.mpl_connect("motion_notify_event", self.handle_mouse_move)

    def handle_key_press(self, event):
        """
        Handles key presses for scaling, clipping, picking modes, etc.
        """
        # ... (Implement key press logic similar to the original Pick function)

    def handle_mouse_click(self, event):
        """
        Handles mouse clicks for adding picks to dictionaries.
        """
        # ... (Implement mouse click logic similar to the original Pick function)

    def handle_mouse_move(self, event):
        """
        Handles mouse movement for highlighting the nearest peak.
        """
def Pick(stream_input, filename="Record Section", outfile="picks.csv"):
    """
    Loads all .su files in a directory into a dictionary of ObsPy Stream objects.
    :param stream_input: The target obspy stream object, composed of all the traces constituting the given shot record
    :param filename: A string representing the name of the target shot record.
    :param outfile: A string representing the filename of the output pick csv file.

    Example usage:
    Pick(ShotRecord, filename="Shot Record", outfile="shotrecord.csv"

    This csv can be imported back into Python as a psPickfile object using pyshot.pick.psPickfile(), or a folder of these csvs may be assimilated into
    a list of pickfile objects using pyshot.pick.assimilate_pickdata()
    """

    # We create global variables both for the stream as we currently are viewing it, and for the stream as it originally
    # was in the data. The first variable 'stream' can be scaled and clipped
    global stream_visual # previously "stream"
    global stream_original
    stream_original = stream_input.copy()
    stream_visual = stream_original.copy()

    # We also create global variables for the filename and directory of the output file
    global outfile_dir
    global outfile_name
    outfile_dir = os.path.dirname(outfile)
    outfile_name = os.path.basename(outfile)


    # First we determine the initial axes of the plot based on the first trace in the Stream.
    # We will use this to determine the x and y limits of the plot.
    mintrace = stream_visual[
        0].stats.su.trace_header.distance_from_center_of_the_source_point_to_the_center_of_the_receiver_group
    maxtrace = stream_visual[
        -1].stats.su.trace_header.distance_from_center_of_the_source_point_to_the_center_of_the_receiver_group
    streamend = np.require(stream_visual[0].stats.endtime, dtype='float32')
    trace_maxvals = []
    for trace in stream_visual:
        trace_maxvals.append(np.max(np.abs(trace.data)))
    maxamp = np.max(trace_maxvals)
    streamfig = plt.figure(figsize=(12, 8))
    stream_visual.plot(type='section', scale=1, fig=streamfig, time_down=True, fillcolors=('blue', 'red'), color='black',
                size=(600, 800))
    ax = plt.gca()
    ax.set_title(filename)
    ax.set_xlabel('Offset (km)')
    ax.set_ylabel('Time (s)')
    ax.set_ylim(streamend, 0)

    global data_delta
    data_delta = stream_original[0].stats.delta

    global scale
    scale = 1

    global clip
    clip = maxamp

    global pickdictPlot1
    pickdictPlot1 = None

    global pickdictPlot2
    pickdictPlot2 = None

    global point
    point, = ax.plot([], [], '_', markersize=25, color='green')

    global pickInfo
    global pickInfo
    pickInfo = ax.text(0.02, 0.02, '', transform=ax.transAxes,
                       verticalalignment='top', bbox=dict(boxstyle='round', facecolor='wheat', alpha=1))

    global pickMode
    pickMode = False

    global pickType
    pickType = 1

    global stdWindow
    stdWindow = 0.02

    offsetlist = []
    for trace in stream_original:
        offsetlist.append(trace.stats.distance / 1000)

    # Precompute every peak in the stream
    peaks = []
    peaktimes = []
    for trace in stream_original:
        peakindices = sp.signal.find_peaks(trace.data)[0]
        troughindices = sp.signal.find_peaks(-trace.data)[0]
        supremaindices = np.concatenate((peakindices, troughindices))
        peaks.append(supremaindices)
        suprematimes_array = trace.stats.delta * supremaindices
        peaktimes.append(suprematimes_array.tolist())

    # Find peak based on cursor coords
    def cursor_peakfinder(event):
        """
        This function finds the peak nearest to the cursor.
        :param event: The cursor event
        :return: The offset, time, amplitude, and standard deviation of the peak nearest to the cursor
        """
        if not event.inaxes:
            return

        nearest_offset = min(offsetlist, key=lambda x: abs(x - event.xdata))
        # This is the offset the cursor is closest to
        # We'll find the index of the trace with this offset
        # and use that to plot the nearest peak
        nearest_index = offsetlist.index(nearest_offset)  # This is the index of the trace the cursor is hovering over
        nearest_trace = stream_original[nearest_index] # This is the trace the cursor is hovering over

        # Find the nearest peak to the cursor
        nearest_peak_time = min(peaktimes[nearest_index], key=lambda y: abs(y - event.ydata))
        global data_delta
        nearest_peak_index = int(nearest_peak_time / data_delta)
        # find the amplitude of the peak nearest to the cursor
        nearest_peak_amplitude = nearest_trace.data[nearest_peak_index]

        # Find the standard deviation window around the peak
        global stdWindow
        n_window = int(stdWindow/2/data_delta)
        window_start = np.max([0, nearest_peak_index-n_window])
        std_window = nearest_trace.data[window_start:nearest_peak_index+n_window]
        std = np.std(std_window)

        return (nearest_offset, nearest_peak_time, nearest_peak_amplitude, std)
    # Highlight nearest peak to the cursor

    def on_move(event):
        """
        This function determines what happens as a function of what the mouse is hovering over
        :param event: The cursor event
        """
        global pickInfo
        if not event.inaxes:
            return
        point.set_data((cursor_peakfinder(event)[0], cursor_peakfinder(event)[1]))
        pickInfo.set_text('Peak amplitude:' + str(cursor_peakfinder(event)[2]) + '\n' + 'Standard deviation:' + str(cursor_peakfinder(event)[3]))
        streamfig.canvas.draw_idle()


    # Add a pick to the pickdict when the user clicks
    pickdict1 = {}
    pickdict2 = {}

    def replot():
        global stream
        old_xlim = ax.get_xlim()
        old_ylim = ax.get_ylim()
        ax.clear()
        stream.plot(type='section', scale=scale, fig=streamfig, time_down=True, fillcolors=('blue', 'red'), color='black',
                    size=(600, 800))
        ax.set_xlim(old_xlim)
        ax.set_ylim(old_ylim)
        ax.set_title(filename)
        global point
        global pickType
        if pickType == 1:
            point, = ax.plot([], [], '_', markersize=20, color='green')
        elif pickType == 2:
            point, = ax.plot([], [], '_', markersize=20, color='fuchsia')
        global pickInfo
        pickInfo = ax.text(0.02, 0.02, '', transform=ax.transAxes,
                           verticalalignment='top', bbox=dict(boxstyle='round', facecolor='wheat', alpha=1))
        global pickdictPlot1
        pickdictPlot1 = ax.scatter(x=list(pickdict1.keys()), y=[pick[0] for pick in pickdict1.values()], marker='X', s=60,
                                   color='white', edgecolors='black')
        global pickdictPlot2
        pickdictPlot2 = ax.scatter(x=list(pickdict2.keys()), y=[pick[0] for pick in pickdict2.values()], marker='X', s=60,
                                   color='black', edgecolors='white')
        streamfig.canvas.draw()

    # Add pick to pickdict on click, and re-plot pickdict every click
    def on_click(event):
        global pickType
        if not event.inaxes:
            return
        pick = cursor_peakfinder(event)
        if pickType == 1:
            # Remove pick if it already exists
            if (pick[0] in pickdict1.keys()):
                if pickdict1[pick[0]] == (pick[1], pick[2], pick[3]):
                    del pickdict1[pick[0]]
                else:
                    add_pick(pick, pickdict1)
            else:
                add_pick(pick, pickdict1)
            global pickdictPlot1
            if pickdictPlot1 is not None:
                pickdictPlot1.remove()
            print(list(pickdict1.keys()))
            print([pick[1] for pick in pickdict1.values()])
        elif pickType == 2:
            # Remove pick if it already exists
            if (pick[0] in pickdict2.keys()):
                if pickdict2[pick[0]] == (pick[1], pick[2], pick[3]):
                    del pickdict2[pick[0]]
                else:
                    add_pick(pick, pickdict2)
            else:
                add_pick(pick, pickdict2)
            global pickdictPlot2
            if pickdictPlot2 is not None:
                pickdictPlot2.remove()
            print(list(pickdict2.keys()))
            print([pick[1] for pick in pickdict2.values()])
        replot()
        streamfig.canvas.draw_idle()
        # ax.plot()

    def on_keypress(event):
        global scale
        global pickMode
        global button
        global motion
        global clip
        global stream
        global stream_original
        global point
        global outfile_dir
        global outfile_name
        if event.key == 'n':
            scale *= 1.1
            replot()
        if event.key == 'm':
            scale /= 1.1
            replot()
        if event.key == 'j':
            clip *= 1.1
            stream = clip_stream(stream_original, clip)
            replot()
        if event.key == 'k':
            clip /= 1.1
            stream = clip_stream(stream_original, clip)
            replot()
        if event.key == 'u':
            clip *= 2
            stream = clip_stream(stream_original, clip)
            replot()
        if event.key == 'i':
            clip /= 2
            stream = clip_stream(stream_original, clip)
            replot()
        if event.key == 'x':
            if pickMode == False:
                pickMode = True
                button = streamfig.canvas.mpl_connect('button_press_event', on_click)
                motion = streamfig.canvas.mpl_connect('motion_notify_event', on_move)
                print('Pick mode on')
            elif pickMode == True:
                pickMode = False
                streamfig.canvas.mpl_disconnect(button)
                streamfig.canvas.mpl_disconnect(motion)
                point, = ax.plot([], [], '_', markersize=25, color='green')
                pickInfo.set_text('')
                replot()
                print('Pick mode off')
        if event.key == '1':
            global pickType
            pickType = 1
            point.set_color('green')
            print('Picking primary arrivals')
        if event.key == '2':
            pickType = 2
            point.set_color('fuchsia')
            print('Picking secondary arrivals')
        if event.key == 'v':
            # Save the picks to a csv file
            output_1 = outfile_dir + '/primary_' + outfile_name
            output_2 = outfile_dir + '/secondary_' + outfile_name
            key_pairs = []
            for key in pickdict1.keys():
                for key2 in pickdict2.keys():
                    if key == key2:
                        key_pairs.append(key)
            with open(output_1, 'w') as f:
                for key in key_pairs:
                    f.write("%s,%s,%s,%s\n" % (key*1000, pickdict1[key][0], pickdict1[key][1], pickdict1[key][2]))
            print('Picks saved to ' + output_1)
            with open(output_2, 'w') as f:
                for key in key_pairs:
                    f.write("%s,%s,%s,%s\n" % (key*1000, pickdict2[key][0], pickdict2[key][1], pickdict2[key][2]))
            print('Picks saved to ' + output_2)
        if event.key == 'c':
            # Save the picks to a csv file
            output_1 = outfile_dir + '/firn_' + outfile_name
            with open(output_1, 'w') as f:
                for key in pickdict1.keys():
                    f.write("%s,%s,%s,%s\n" % (key*1000, pickdict1[key][0], pickdict1[key][1], pickdict1[key][2]))
            print('Picks saved to ' + output_1)
    streamfig.canvas.mpl_connect('key_press_event', on_keypress)
    plt.show(block=True)