import os
import obspy as op


def loadshots(directory_path):
    """
    Loads all .su files in a directory into a dictionary of ObsPy Stream objects.
    """
    streamdict = {}  # Initialize an empty dictionary.
    # Loop over all files in the directory.
    for filename in os.listdir(directory_path):
        if filename.endswith(".su"):  # Make sure we're only reading .su files.

            # Generate the full file path.
            full_path = os.path.join(directory_path, filename)

            # Read the .su file into an ObsPy Stream object.
            stream = op.read(full_path, format="SU", unpack_trace_headers=True)

            # Use SU headers to set the ObsPy trace headers
            for trace in stream:
                # trace.stats.distance = abs(trace.stats.su.trace_header.distance_from_center_of_the_source_point_to_the_center_of_the_receiver_group)
                trace.stats.distance = trace.stats.su.trace_header.distance_from_center_of_the_source_point_to_the_center_of_the_receiver_group
                # if trace.stats.distance == 0:
                #     stream.remove(trace)
            # Store the Stream object in the dictionary.
            streamdict[filename] = stream

    return streamdict
