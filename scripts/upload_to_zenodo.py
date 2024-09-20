# Description: This script is used to upload a file to an existing Zenodo deposition.
import requests
import os
import sys


# this works!

access_token = '<access-token>'
headers = {'Authorization': f'Bearer {access_token}'}

deposition_id = '<deposition-id>'

url = f'https://zenodo.org/api/deposit/depositions/{deposition_id}'
response = requests.get(url, headers=headers)

deposition_data = response.json()
bucket_url = deposition_data['links']['bucket']

paths = [
    # '/work/FAC/FBM/DMF/pengel/spirit/aprasad/BACKUP_current/20230313_apis_species_comparison_zips/README.md',
    # '/work/FAC/FBM/DMF/pengel/spirit/aprasad/BACKUP_current/20230313_apis_species_comparison_zips/05_assembly.zip',
    # '/work/FAC/FBM/DMF/pengel/spirit/aprasad/BACKUP_current/20230313_apis_species_comparison_zips/06_metagenomicORFs.zip',
    # '/work/FAC/FBM/DMF/pengel/spirit/aprasad/BACKUP_current/20230313_apis_species_comparison_zips/11_phylogenies.zip',
    '/work/FAC/FBM/DMF/pengel/spirit/aprasad/BACKUP_current/20230313_apis_species_comparison_zips/figures.zip',
    '/work/FAC/FBM/DMF/pengel/spirit/aprasad/BACKUP_current/20230313_apis_species_comparison_zips/09_MAGs_collection.zip'
    ]


class ProgressFile:
    def __init__(self, file_path):
        self.file_path = file_path
        self.file = open(file_path, 'rb')
        self.total_size = os.path.getsize(file_path)  # Get the size of the file
        self.uploaded_size = 0
        
    def __len__(self):
        return self.total_size
        
    def read(self, block_size):
        data = self.file.read(block_size)
        if data:
            self.uploaded_size += len(data)
            percent_complete = (self.uploaded_size / self.total_size) * 100
            # Print progress every 5%
            if int(percent_complete) % 5 == 0:
                sys.stdout.write(f'\rUpload progress for {os.path.basename(self.file_path)}: {percent_complete:.2f}%')
                sys.stdout.flush()
        return data

for path in paths:
    progress_file = ProgressFile(path)
    upload_url = f'{bucket_url}/{os.path.basename(path)}'
    upload_response = requests.put(upload_url, headers=headers, data=progress_file)
    progress_file.file.close()  # Close the file after uploading
    break
    
    

# url = f'https://zenodo.org/api/deposit/depositions/{deposition_id}'
# response = requests.get(url, headers=headers)

# # Check if the deposition details were retrieved successfully
# if response.status_code != 200:
#     print(f"Error retrieving deposition: {response.status_code}")
#     print(response.json())
# else:
#     # Parse the response JSON
#     deposition_data = response.json()

#     # Check if there are any files associated with the deposition
#     if 'files' in deposition_data and deposition_data['files']:
#         print("List of files in the deposition:")
#         for file_info in deposition_data['files']:
#             file_name = file_info['filename']
#             file_size = file_info['filesize']
#             print(f"File: {file_name}, Size: {file_size} bytes")
#     else:
#         print("No files have been uploaded to this deposition.")