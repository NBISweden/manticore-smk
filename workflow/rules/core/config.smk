def get_app_params(wildcards, rule):
    return config

def read_pairs_dataframe():
    """List reads as paired-end or single-end"""
    df = reads.groupby(level=["SM", "unit"]).size().to_frame("pe")
    df = reads[reads["id"] == 1].droplevel("id").join(df["pe"])
    df["pe"] = df["pe"].map({1: "se", 2: "pe"})
    return df
