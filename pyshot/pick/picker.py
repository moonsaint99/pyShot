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

class Picker:
    def __init__(self, stream_input, filename="Record Section", outfile="picks.csv"):
        """
        Initializes the Picker class for interactive seismic data picking.
        
        :param stream_input: The target ObsPy Stream object, composed of all the traces constituting the given shot record
        :param filename: A string representing the name of the target shot record.
        :param outfile: A string representing the filename of the output pick csv file.
        """
        self.stream_original = stream_input.copy()
        self.stream_visual = self.stream_original.copy()
        self.filename = filename
        self.outfile_dir = os.path.dirname(outfile)
        self.outfile_name = os.path.basename(outfile)
        
        self.data_delta = self.stream_original[0].stats.delta
        self.scale = 1
        self.clip = self._calculate_max_amplitude()
        
        self.pickdict1 = {}
        self.pickdict2 = {}
        self.pickMode = False
        self.pickType = 1
        self.stdWindow = 0.02
        
        self._initialize_plot()
        self._precompute_peaks()
        
    def _calculate_max_amplitude(self):
        return max(np.max(np.abs(trace.data)) for trace in self.stream_visual)
    
    def _initialize_plot(self):
        self.fig, self.ax = plt.subplots(figsize=(12, 8))
        self.stream_visual.plot(type='section', scale=self.scale, fig=self.fig, time_down=True,
                                fillcolors=('blue', 'red'), color='black', size=(600, 800))
        self._set_plot_properties()
        
        self.point, = self.ax.plot([], [], '_', markersize=25, color='green')
        self.pickInfo = self.ax.text(0.02, 0.02, '', transform=self.ax.transAxes,
                                     verticalalignment='top', bbox=dict(boxstyle='round', facecolor='wheat', alpha=1))
        
        self.fig.canvas.mpl_connect('key_press_event', self.on_keypress)
        self.fig.canvas.mpl_connect('button_press_event', self.on_click)
        self.fig.canvas.mpl_connect('motion_notify_event', self.on_move)
    
    def _set_plot_properties(self):
        self.ax.set_title(self.filename)
        self.ax.set_xlabel('Offset (km)')
        self.ax.set_ylabel('Time (s)')
        self.ax.set_ylim(np.require(self.stream_visual[0].stats.endtime, dtype='float32'), 0)
    
    def _precompute_peaks(self):
        self.offsetlist = [trace.stats.distance / 1000 for trace in self.stream_original]
        self.peaks = []
        self.peaktimes = []
        for trace in self.stream_original:
            peakindices = sp.signal.find_peaks(trace.data)[0]
            troughindices = sp.signal.find_peaks(-trace.data)[0]
            supremaindices = np.concatenate((peakindices, troughindices))
            self.peaks.append(supremaindices)
            self.peaktimes.append((trace.stats.delta * supremaindices).tolist())
    def cursor_peakfinder(self, event):
        if not event.inaxes or not self.pickMode:
            return None
        
        nearest_offset = min(self.offsetlist, key=lambda x: abs(x - event.xdata))
        nearest_index = self.offsetlist.index(nearest_offset)
        nearest_trace = self.stream_original[nearest_index]
        
        nearest_peak_time = min(self.peaktimes[nearest_index], key=lambda y: abs(y - event.ydata))
        nearest_peak_index = self.peaks[nearest_index][self.peaktimes[nearest_index].index(nearest_peak_time)]
        nearest_peak_amplitude = nearest_trace.data[nearest_peak_index]
        
        n_window = int(self.stdWindow / self.data_delta)
        window_start = max(0, nearest_peak_index - n_window // 2)
        window_end = min(len(nearest_trace.data), nearest_peak_index + n_window // 2)
        std_window = nearest_trace.data[window_start:window_end]
        std = np.std(std_window)
        
        return (nearest_offset, nearest_peak_time, nearest_peak_amplitude, std)
    
    def on_move(self, event):
        if not self.pickMode:
            self.point.set_data([], [])
            self.pickInfo.set_text('')
            self.fig.canvas.draw_idle()
            return

        peak_info = self.cursor_peakfinder(event)
        if peak_info:
            self.point.set_data([peak_info[0]], [peak_info[1]])  # Wrap single values in lists
            self.pickInfo.set_text(f'Peak amplitude: {peak_info[2]:.2f}\nStandard deviation: {peak_info[3]:.2f}')
            self.fig.canvas.draw_idle()
    
    def on_click(self, event):
        if not self.pickMode or not event.inaxes:
            return
        
        pick = self.cursor_peakfinder(event)
        if not pick:
            return
        
        current_dict = self.pickdict1 if self.pickType == 1 else self.pickdict2
        
        if pick[0] in current_dict and current_dict[pick[0]] == pick[1:]:
            del current_dict[pick[0]]
        else:
            current_dict[pick[0]] = pick[1:]
        
        self.replot()
    
    def replot(self):
        # Store current xlim, ylim, and title
        xlim = self.ax.get_xlim()
        ylim = self.ax.get_ylim()
        title = self.ax.get_title()

        self.ax.clear()
        self.stream_visual.plot(type='section', scale=self.scale, fig=self.fig, time_down=True,
                                fillcolors=('blue', 'red'), color='black', size=(600, 800))
        self._set_plot_properties()
        
        self.point, = self.ax.plot([], [], '_', markersize=20, color='green' if self.pickType == 1 else 'fuchsia')
        self.pickInfo = self.ax.text(0.02, 0.02, '', transform=self.ax.transAxes,
                                     verticalalignment='top', bbox=dict(boxstyle='round', facecolor='wheat', alpha=1))
        
        self.ax.scatter(x=list(self.pickdict1.keys()), y=[pick[0] for pick in self.pickdict1.values()],
                        marker='X', s=60, color='white', edgecolors='black')
        self.ax.scatter(x=list(self.pickdict2.keys()), y=[pick[0] for pick in self.pickdict2.values()],
                        marker='X', s=60, color='black', edgecolors='white')
        
        # Restore previous xlim, ylim, and title
        self.ax.set_xlim(xlim)
        self.ax.set_ylim(ylim)
        self.ax.set_title(title)
        
        self.fig.canvas.draw()
    
    def on_keypress(self, event):
        key_actions = {
            'n': lambda: self._adjust_scale(1.1),
            'm': lambda: self._adjust_scale(1/1.1),
            'j': lambda: self._adjust_clip(1.1),
            'k': lambda: self._adjust_clip(1/1.1),
            'u': lambda: self._adjust_clip(2),
            'i': lambda: self._adjust_clip(1/2),
            'x': self._toggle_pick_mode,
            '1': lambda: self._set_pick_type(1),
            '2': lambda: self._set_pick_type(2),
            'v': self._save_picks,
            'c': self._save_firn_picks
        }
        
        if event.key in key_actions:
            key_actions[event.key]()
    
    def _adjust_scale(self, factor):
        self.scale *= factor
        self.replot()
    
    def _adjust_clip(self, factor):
        self.clip *= factor
        self.stream_visual = self._clip_stream(self.stream_original, self.clip)
        self.replot()
    
    def _clip_stream(self, stream, clip_value):
        clipped_stream = stream.copy()
        for trace in clipped_stream:
            trace.data = np.clip(trace.data, -clip_value, clip_value)
        return clipped_stream
    
    def _toggle_pick_mode(self):
        self.pickMode = not self.pickMode
        print(f"Pick mode {'on' if self.pickMode else 'off'}")
        if not self.pickMode:
            self.point.set_data([], [])
            self.pickInfo.set_text('')
            self.replot()
    
    def _set_pick_type(self, pick_type):
        self.pickType = pick_type
        self.point.set_color('green' if pick_type == 1 else 'fuchsia')
        print(f"Picking {'primary' if pick_type == 1 else 'secondary'} arrivals")
    
    def _save_picks(self):
        self._save_pick_file('primary', self.pickdict1)
        self._save_pick_file('secondary', self.pickdict2)
    
    def _save_firn_picks(self):
        self._save_pick_file('firn', self.pickdict1)
    
    def _save_pick_file(self, prefix, pickdict):
        output_file = os.path.join(self.outfile_dir, f'{prefix}_{self.outfile_name}')
        with open(output_file, 'w') as f:
            for key, value in pickdict.items():
                f.write(f"{key*1000},{value[0]},{value[1]},{value[2]}\n")
        print(f'Picks saved to {output_file}')
    
    def run(self):
        plt.show(block=True)

def Pick(stream_input, filename="Record Section", outfile="picks.csv"):
    """
    Creates and runs a Picker instance for interactive seismic data picking.
    
    :param stream_input: The target ObsPy Stream object, composed of all the traces constituting the given shot record
    :param filename: A string representing the name of the target shot record.
    :param outfile: A string representing the filename of the output pick csv file.
    
    Example usage:
    Pick(ShotRecord, filename="Shot Record", outfile="shotrecord.csv")
    
    This csv can be imported back into Python as a psPickfile object using pyshot.pick.psPickfile(),
    or a folder of these csvs may be assimilated into a list of pickfile objects using pyshot.pick.assimilate_pickdata()
    """
    picker = Picker(stream_input, filename, outfile)
    picker.run()