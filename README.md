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

4. Use a command like `aws s3 sync path/to/your/dataset/ s3://veda-data-store-staging/EIS/<your-dataset-name>/` to push the data into S3.
  Note the location in S3 (`EIS` "subdirectory").
  Note also the trailing slashes for the sync command -- this will sync the _contents_ of `path/to/your/dataset/` with the contents of `.../EIS/<your-dataset-name>/`
  (i.e., it will put the file `path/to/your/dataset/somefile` into `.../EIS/<your-dataset-name>/somefile`.
  If you omit the trailing slashes, you will instead get something like `.../EIS/<your-dataset-name>/dataset/somefile`.
  ).

  <details>

  The following slightly simpler version of this code can be used to get S3 credentials for any interactive environment that uses IAM-based authentication (including MAAP and the EIS SMCE DaskHubs):

  ```python
  import boto3
  session = boto3.Session()
  cred = session.get_credentials().get_frozen_credentials()
  print(f"""
  export AWS_ACCESS_KEY_ID={cred.access_key}
  export AWS_SECRET_ACCESS_KEY={cred.secret_key}
  export AWS_SESSION_TOKEN={cred.token}
  """)
  ```

  Similar to the instructions above, this will generate environment variables that can be transferred to any machine to give that machine temporary permissions equivalent to the AWS system on which that code was originally run (i.e., a machine with these variables set literally "pretends" to be the same as the original machine, from AWS's perspective).
  These permissions can also be used to download datasets; e.g., `aws s3 cp s3://veda-data-store-staging/path/to/your/file ./local/destination/path/to/file` or `aws s3 sync s3://veda-data-store-staging/path/to/some/folder/ ./local/folder` (note, first argument is the source, the second is the destination).
  However, **please use this download capability carefully and sparingly**, as moving data out of AWS incurs egress costs, which are low for small volumes of data but can quickly add up for large volumes.
  For tests on small/medium-sized datasets, it's totally fine, but please avoid downloading large datasets in this way.

  </details>

## Processed datasets

- [SPL3SMP](SPL3SMP/README.md) -- SMAP L3 Radiometer Global Daily 36 km EASE-Grid Soil moisture
- [FIRMS active fires](FIRMS/README.md) -- Active fire data from MODIS (MOD14) and VIIRS (VNP14) from FIRMS/LANCE.
