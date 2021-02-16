# JERICHO-E-usage data package

The JERICHO-E-usage dataset including code for compilation. The dataset includes time series of useful energy consumption patterns for energy system modeling. The time series data covers 38 NUTS2 regions in germany and has an hourly resolution. The dataset distinguished between four sectors - the residential, the industrial, the commerce, and the mobility sectors. Useful energy types comprise space heating, warm water, process heating, space cooling, process cooling, mechanical energy, information and communication technology, and lighting.

See <Link to publication> for a detailed description of the data and the underlying methodology.

## Usage Instructions

### Install environment

The most convenient way to use this notebook is to use a conda environment: 
<ol>
    <li>Install conda: https://docs.conda.io/projects/conda/en/latest/user-guide/install/</li>
    <li>Download the Git repository</li>
    <li>Create the environment using the environment.yml file provided in the Git reposity by executing the following command in the Git repository:
        <ul>
        <li><code>conda env create -f environment.yml</code></li>
        </ul>
    </li>
    <li>Activate the environment by executing the following command:
        <ul>
        <li><code>conda activate JERICHO-E-usage</code></li>
        </ul>
    </li>
    <li>Start Jupyter Notebook by executing the following command in the Git repository
        <ul>
        <li><code>jupyter notebook</code></li>
        </ul>
    </li>
</ol> 

### Run scripts

For easy use of the scripts we have added a Jupyter Notebook with further instructions on the workflow.

## License

This repository is published under the [MIT License](LICENSE.md).
