install.packages("base64enc")

library(base64enc)

# Read the .dll file and encode it to Base64
dll_path <- "/Users/youbeenshim/Downloads/trajv9-64"
output_path <- ""
getwd()
# Encode the file
base64encode(dll_path, output_path)

library(base64enc)

# Read the encoded text file and decode back to .dll
encoded_path <- "/Users/youbeenshim/Downloads/trajv9-64/traj.dll"
output_dll <- ""

# Decode the file
base64decode(file = encoded_path, output = output_dll)



# Read the binary file
dll_data <- readBin("C:/path/to/traj.dll", "raw", n = file.info("C:/path/to/traj.dll")$size)

# Encode to Base64
encoded_string <- base64enc::base64encode(dll_data)

# Write to text file
writeLines(encoded_string, "C:/path/to/traj_encoded.txt")
# Read the encoded text
encoded_string <- readLines("C:/path/to/traj_encoded.txt", warn = FALSE)
encoded_string <- paste(encoded_string, collapse = "")

# Decode from Base64
dll_data <- base64enc::base64decode(encoded_string)

# Write binary file
writeBin(dll_data, "C:/path/to/traj.dll")

