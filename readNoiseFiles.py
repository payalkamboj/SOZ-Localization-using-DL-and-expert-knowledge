#Automated code for Noise list from DL part of DeepXSOZ
import re

class extractNoiseLabels():
    def __init__(self, folder):
         self.filename = folder
         self.extractNoiseICNumbers()

    def extractNoiseICNumbers(self):

        filename = self.filename  # Replace with the actual file path and name

        try:
            with open(filename, 'r') as file:
                lines = file.readlines()
                # Filter lines that start with "Noise"
                filtered_lines = [line.strip() for line in lines if line.startswith('Noise')]
                extracted_numbers = []

                for line in filtered_lines:
                    parts = line.split('_')
                    # Find the part that starts with "NoiseASUAI"
                    noise_part = next((part for part in parts if part.startswith('0')), None)
                    if noise_part:
                         try:
                              number = int(noise_part[len('0'):])
                              extracted_numbers.append(number)
                         except ValueError:
                             pass
  

                print(extracted_numbers)
                unique_numbers_set = set(extracted_numbers)
                print(unique_numbers_set)

                for number in unique_numbers_set:
                    list = []
                    print(f"Lines for number {number}:")
                    for line in filtered_lines:
             
                      if(f"_00{number}_" in line or f"_0{number}_" in line) and " 0.0" in line:
                        print(line)
                        pattern = r'_([0-9]+)_\.png'
                        match = re.search(pattern, line)
                        if match:
                            ic_Number = match.group(1)
                            list.append(ic_Number)
                            print(match.group(1))
                    # Convert the list to a string with '[' and ']' included
                    list_str = '[' + '; '.join(map(str, list)) + ']'

                    file_name = f"Noise{number}.txt"
                    with open(file_name, 'w') as file:
                         file.write(list_str)



        except FileNotFoundError:
            print(f"File '{filename}' not found.")
        except PermissionError:
            print(f"Permission denied for file '{filename}'.")
        except Exception as e:
            print(f"An error occurred: {str(e)}")
