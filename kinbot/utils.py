import time


def tail(file_path, lines=10):
    """Returns the specified last number of lines of a file.

    @param file_path: The n to retrieve the last lines from.
    @param lines: The number of lines to be retrieved.
    @return str: The last number of lines
    """
    with open(file_path, 'br') as f:
        total_lines_wanted = lines

        block_size = 1024
        f.seek(0, 2)
        block_end_byte = f.tell()
        lines_to_go = total_lines_wanted
        block_number = -1
        blocks = []
        while lines_to_go > 0 and block_end_byte > 0:
            if block_end_byte - block_size > 0:
                f.seek(block_number * block_size, 2)
                blocks.append(f.read(block_size))
            else:
                f.seek(0, 0)
                blocks.append(f.read(block_end_byte))
            lines_found = blocks[-1].count(b'\n')
            lines_to_go -= lines_found
            block_end_byte -= block_size
            block_number -= 1
        all_read_text = b''.join(reversed(blocks))
    return b'\n'.join(all_read_text.splitlines()[-total_lines_wanted:]).decode()


def iowait(logfile, qc_code):
    """Wait for I/O to be done to a specific file.
    There is a maxtime.
    """
    maxtime = 10  # s
    clock = 0

    if qc_code == 'gauss':
        endings = ['Normal termination', 'Error termination']
    elif qc_code == 'qchem':
        endings = ['Thank you very much for using Q-Chem', 'Please submit a crash report']
    else:
        endings = []
    while True:
        logtail = tail(logfile)
        if any([ee in logtail for ee in endings]):
            break
        clock += 1
        if clock > maxtime:
            break
        time.sleep(1)
