import argparse
import os
import pandas as pd
import logging


def setup_logging(log_file):
    """Configure logging settings."""
    logging.basicConfig(filename=log_file,
                        level=logging.INFO,
                        format='%(asctime)s - %(levelname)s - %(message)s')
    # If you also want to print the log messages to the console, uncomment the following lines:
    # console_handler = logging.StreamHandler()
    # logging.getLogger().addHandler(console_handler)

def validate_table(table: pd.DataFrame, table_name: str, CDE: pd.DataFrame):

    retval = 1

    # Filter out rows specific to the given table_name from the CDE
    specific_cde_df = CDE[CDE['Table'] == table_name]
    
    # Extract fields that have a data type of "Enum" and retrieve their validation entries
    enum_fields_dict = dict(zip(specific_cde_df[specific_cde_df['DataType'] == "Enum"]['Field'], 
                               specific_cde_df[specific_cde_df['DataType'] == "Enum"]['Validation']))
    
    # Extract fields that are marked as "Required"
    required_fields = specific_cde_df[specific_cde_df['Required'] == "Required"]['Field'].tolist()
    optional_fields = specific_cde_df[specific_cde_df['Required'] == "Optional"]['Field'].tolist()

    table = force_enum_string(table, table_name, CDE)

    # Check for missing "Required" fields
    missing_required_fields = [field for field in required_fields if field not in table.columns]
    
    if missing_required_fields:
        logging.info(f"\tMissing Required Fields in {table_name}: {', '.join(missing_required_fields)}")
    else:
        logging.info(f"\tAll required fields are present in {table_name}.")

    # Check for empty or NaN values
    for test_field,test_name in zip([required_fields, optional_fields], ["Required", "Optional"]):
        empty_or_nan_fields = {}
        for field in test_field:
            if field in table.columns:
                invalid_count = table[field].isna().sum()
                if invalid_count > 0:
                    empty_or_nan_fields[field] = invalid_count
                    
        if empty_or_nan_fields:
            logging.info(f"\t\t{test_name} Fields with Empty or NaN values:")
            for field, count in empty_or_nan_fields.items():
                logging.info(f"\t\t\t- {field}: {count} rows")
            retval = 0
        else:
            logging.info(f"No empty or NaN values found in {test_name} fields.")
    


    # Check for invalid Enum field values
    invalid_field_values = {}
    for field, validation_str in enum_fields_dict.items():
        valid_values = eval(validation_str)
        if field in table.columns:
            invalid_values = table[~table[field].isin(valid_values)][field].unique()
            if invalid_values.any():
                invalid_field_values[field] = invalid_values
    
    if invalid_field_values:
        logging.info("\tInvalid Field/Value pairs:")
        for field, values in invalid_field_values.items():
            logging.info(f"\t\t\t- {field}: {', '.join(map(str, values))}")
        retval = 0
    else:
        logging.info(f"\tAll Enum fields have valid values in {table_name}.")

    return retval

######## HELPERS ########
# Define a function to only capitalize the first letter of a string
def capitalize_first_letter(s):
    if not isinstance(s, str) or len(s) == 0:  # Check if the value is a string and non-empty
        return s
    return s[0].upper() + s[1:]

def force_enum_string(df, df_name, CDE):

    string_enum_fields = CDE[(CDE["Table"] == df_name) & 
                                (CDE["DataType"].isin(["Enum", "String"]))]["Field"].tolist()
    # Convert the specified columns to string data type using astype() without a loop
    columns_to_convert = {col: 'str' for col in string_enum_fields if col in df.columns}
    df = df.astype(columns_to_convert)

    # enum_fields = CDE[ (CDE["Table"] == df_name) & 
    #                             (CDE["DataType"]=="Enum") ]["Field"].tolist()
    
    for col in string_enum_fields:
        if col in df.columns and col not in ["assay", "file_type"]:
            df[col] = df[col].apply(capitalize_first_letter)

    return df

def validate_file(table_name, path):
    # Construct the path to CSD.csv
    cde_file_path = os.path.join(path, "CDE.csv")
    ref_file_path = os.path.join(path, f"{table_name}.csv")
    logging.info(f"## {table_name} ({ref_file_path})")

    # Load the CDE.csv file and the reference table
    CDE_df = pd.read_csv(cde_file_path)

    reference_df = pd.read_csv(ref_file_path)

    retval = validate_table(reference_df, table_name, CDE_df)
    return retval




def main():
    # Set up the argument parser
    parser = argparse.ArgumentParser(description="A command-line tool to validate against CSD.csv.")
    
    # Add arguments
    parser.add_argument("table", choices=["SAMPLE", "STUDY", "PROTOCOL", "CLINPATH", "SUBJECT", "CDE", "ALL"],
                        help="Specify the table for validation.")
    parser.add_argument("--path", default=os.getcwd(),
                        help="Path to the directory containing CSD.csv. Defaults to the current working directory.")
    parser.add_argument("--logfile", default="validation.log",
                        help="Name of the logfile. Defaults to 'validation.log'.")
    
    # Parse the arguments
    args = parser.parse_args()

    # Set up logging to a file named 'validation.log' (or choose another name)
    setup_logging(args.logfile)

    # Call the validation function
    if args.table == "ALL":
        for table_name in ["SAMPLE", "STUDY", "PROTOCOL", "CLINPATH", "SUBJECT"]:
            retval = validate_file(table_name, args.path)
    else:
        retval = validate_file(args.table, args.path)
    

if __name__ == "__main__":
    main()
