# veda-data-processing

Scripts for data downloading, transformation, and related processing for the VEDA project.

## Instructions for uploading data

1. On the machine where your data are, make sure you have the AWS CLI installed.
  Note that it's a pre-compiled binary, so you don't need admin permissions to do this (but you do need to make sure wherever you put the `aws` binary is on your `PATH`).

2. Obtain AWS credentials from MAAP.
  In a Jupyter notebook running in MAAP, run the following code:

  ```python
  import boto3
  sts_client = boto3.client('sts')
  assumed_role_object = sts_client.assume_role(
      RoleArn='arn:aws:iam::114506680961:role/MAAP-users-VEDA-role',
      RoleSessionName='TestSession'
  )
  credentials = assumed_role_object['Credentials']
  print(f"""
  # Credentials expire: {credentials["Expiration"].strftime("%Y-%m-%d %H:%M")}
  export AWS_ACCESS_KEY_ID={credentials["AccessKeyId"]}
  export AWS_SECRET_ACCESS_KEY={credentials["SecretAccessKey"]}
  export AWS_SESSION_TOKEN={credentials["SessionToken"]}
  """)
  ```
    
3. Open a terminal on the machine where the data are stored.
  Copy the output of that command, paste it into that terminal, and hit Enter.
  This will set your temporary AWS credentials.
  To confirm that the variables were set, you can run `env | grep AWS`.
  To confirm that you can actually access the relevant bucket, run `aws s3 ls s3://veda-data-store-staging/`.

4. Use a command like `aws s3 sync s3://veda-data-store-staging/EIS/<your-dataset-name>/ path/to/your/dataset/` to push the data into S3.
  Note the location in S3 (`EIS` "subdirectory").
  Note also the trailing slashes for the sync command -- this will sync the _contents_ of `path/to/your/dataset/` with the contents of `.../EIS/<your-dataset-name>/`
  (i.e., it will put the file `path/to/your/dataset/somefile` into `.../EIS/<your-dataset-name>/somefile`.
  If you omit the trailing slashes, you will instead get something like `.../EIS/<your-dataset-name>/dataset/somefile`.
  ).

## Processed datasets

- [SPL3SMP](SPL3SMP/README.md) -- SMAP L3 Radiometer Global Daily 36 km EASE-Grid Soil moisture
- [FIRMS active fires](FIRMS/README.md) -- Active fire data from MODIS (MOD14) and VIIRS (VNP14) from FIRMS/LANCE.
