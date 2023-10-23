The file checksum_sha256.txt is included in this zip package to be able to verify that the csv file is
complete and doesn't contain any errors after the download. It includes a SHA256 checksum which is a sequence of numbers and
letters and acts as a fingerprint for the csv file.

How to verify on Windows 10
1. Extract the csv file from the zip package. Detailed instructions:
        Right-click the downloaded zip package and choose Extract all. You can select any destination folder.
        Open the destination folder after extracting.

2. Open a command window in the extracted folder. Detailed instructions:
        Ensure the extracted folder is open in the File Explorer, and that the File Explorer window is active.
        Click the address field of the File Explorer window (the path to the current folder should become highlighted).
        Type "cmd" without enclosing quotes and press Enter.
        A command windows should now be opened.

3.Type "CertUtil -hashfile" without enclosing quotes in the terminal window and then press the space button

4.Type the csv file name, or type only the first few letters in the file name and press the tab key. The file name will be filled 
   in automatically when you press tab.

5.Press the space button, type "SHA256" without quotes and press Enter

6.Compare checksums.Detailed instructions:
        The SHA256 checksum that now is shown in the command window should be the same as the checksum in checksum_sha256.txt


How to verify on macOS
1.Extract the csv file from the zip package

2.Open a terminal window.Detailed instructions:
        Use keyboard combination "Command + Space" to open Spotlight.Write "Terminal" without
        enclosing quotes in the Spotlight search field and press Enter

3.Type "shasum -a 256" without enclosing quotes in the terminal window and then press the space button

4.Drag the csv file into the terminal window. This will generate the file path in the terminal window so you don't have to type it.

5.Press Enter

6.Compare checksums:
        The SHA256 checksum that now is shown in the terminal should be the same as the checksum in checksum_sha256.txt

How to verify on Linux
1.Open a terminal

2.Enter the following command:

        sha256sum < path to csv file >

3.Press Enter

4.The SHA256 checksum that now is shown in the terminal should be the same as the checksum in checksum_sha256.txt
