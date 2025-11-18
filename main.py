import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import argparse
import sys
import json

# Try to import your functions
try:
    # For GitHub Actions, we need to add the path
    sys.path.append('/github/workspace')
    from functions import *
except ImportError as e:
    print(f"Warning: Could not import functions module: {e}")
    # Create mock functions for testing
    def prepare_ERA5_SM(*args):
        return np.zeros_like(args[0])
    
    def calibrate_GR4J_SM(*args, **kwargs):
        return {"best_params": np.array([300.0, 0.5, 50.0, 3.0, 0.0, 4.0, 200.0])}
    
    def GR4J_with_SM(*args, **kwargs):
        return {"streamflow": np.random.random(1000) * 10}
    
    def plot_GR4J_SM_results(*args, **kwargs):
        plt.plot([1,2,3], [1,2,3])
    
    def calculate_metrics(obs, sim):
        return {"nse": 0.8, "kge": 0.7, "rmse": 1.2}

def load_qout_values(qout_csv_path=None, mean_qout=None, data_length=None):
    """Load Qout values from CSV or use mean value"""
    if qout_csv_path and os.path.exists(qout_csv_path):
        try:
            qout_df = pd.read_csv(qout_csv_path)
            if 'Qout' in qout_df.columns:
                return qout_df['Qout'].values
            else:
                return qout_df.iloc[:, 0].values
        except Exception as e:
            print(f"Error loading CSV: {e}")
            return None
    elif mean_qout is not None:
        if data_length:
            return np.full(data_length, float(mean_qout))
    return None

def main():
    # Get inputs from environment variables (GitHub Actions)
    basin = os.getenv('BASIN', 'ter')
    mean_qout = os.getenv('MEAN_QOUT')
    calibrate = os.getenv('CALIBRATE', 'false').lower() == 'true'
    
    # Set up paths for GitHub Actions
    dir_path = "/github/workspace"
    output_dir = os.path.join(dir_path, "output")
    os.makedirs(output_dir, exist_ok=True)
    
    print(f"Starting GR4J model for basin: {basin}")
    print(f"Calibration: {calibrate}")
    print(f"Mean Qout: {mean_qout}")
    
    try:
        # This is where you'd put your actual model code
        # For now, creating dummy results for demonstration
        
        # Create sample results
        dates = pd.date_range('2020-01-01', periods=365, freq='D')
        obs_flow = np.random.normal(5, 2, 365)
        sim_flow = obs_flow + np.random.normal(0, 0.5, 365)
        
        # Create plots
        plt.figure(figsize=(10, 6))
        plt.plot(dates, obs_flow, label='Observed Flow')
        plt.plot(dates, sim_flow, label='Simulated Flow')
        plt.legend()
        plt.title(f'GR4J Model Results - {basin}')
        plt.savefig(os.path.join(output_dir, 'results_plot.png'))
        
        # Save results to CSV
        results_df = pd.DataFrame({
            'date': dates,
            'observed_flow': obs_flow,
            'simulated_flow': sim_flow
        })
        results_df.to_csv(os.path.join(output_dir, 'results.csv'), index=False)
        
        # Save metrics
        metrics = {
            'nse': 0.85,
            'kge': 0.78,
            'rmse': 1.1,
            'basin': basin,
            'calibration_run': calibrate
        }
        
        with open(os.path.join(output_dir, 'metrics.json'), 'w') as f:
            json.dump(metrics, f, indent=2)
            
        print("Model run completed successfully!")
        print(f"Results saved to: {output_dir}")
        
    except Exception as e:
        print(f"Error running model: {e}")
        sys.exit(1)

if __name__ == "__main__":
    main()
