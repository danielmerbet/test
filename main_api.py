from flask import Flask, request, jsonify, send_file
from flask_cors import CORS
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import json
import os
import threading
import time
from datetime import datetime, timedelta

app = Flask(__name__)
CORS(app)

# Global variable to store simulation status and results
simulation_status = {
    'running': False,
    'progress': 0,
    'message': '',
    'results': None
}

def run_gr4j_simulation(params):
    """Run the GR4J simulation with progress updates"""
    global simulation_status
    
    try:
        simulation_status['running'] = True
        simulation_status['progress'] = 0
        simulation_status['message'] = 'Starting simulation...'
        
        # Simulate your actual GR4J model here
        time.sleep(2)
        simulation_status['progress'] = 25
        simulation_status['message'] = 'Loading input data...'
        
        # Generate sample data
        dates = pd.date_range('2020-01-01', periods=365, freq='D')
        obs_flow = np.random.normal(5, 2, 365)
        sim_flow = obs_flow + np.random.normal(0, 0.5, 365)
        
        time.sleep(2)
        simulation_status['progress'] = 50
        simulation_status['message'] = 'Running hydrological model...'
        
        # Calculate metrics
        nse = 1 - np.sum((obs_flow - sim_flow)**2) / np.sum((obs_flow - np.mean(obs_flow))**2)
        kge = 0.78  # Simplified
        rmse = np.sqrt(np.mean((obs_flow - sim_flow)**2))
        pbias = np.sum(sim_flow - obs_flow) / np.sum(obs_flow) * 100
        
        time.sleep(2)
        simulation_status['progress'] = 75
        simulation_status['message'] = 'Generating plots and results...'
        
        # Generate volume data
        volume = 1e6 + np.cumsum(sim_flow * 1000)  # Simplified volume calculation
        
        # Create plots
        plt.figure(figsize=(10, 6))
        plt.plot(dates, obs_flow, label='Observed Flow')
        plt.plot(dates, sim_flow, label='Simulated Flow')
        plt.legend()
        plt.title('GR4J Model Results')
        plt.savefig('output/results_plot.png')
        plt.close()
        
        # Prepare results
        results = {
            'metrics': {
                'nse': float(nse),
                'kge': float(kge),
                'rmse': float(rmse),
                'pbias': float(pbias)
            },
            'dates': [d.strftime('%Y-%m-%d') for d in dates],
            'observed_flow': obs_flow.tolist(),
            'simulated_flow': sim_flow.tolist(),
            'simulated_volume': volume.tolist(),
            'parameters': params
        }
        
        simulation_status['progress'] = 100
        simulation_status['message'] = 'Simulation completed successfully!'
        simulation_status['results'] = results
        
        # Save results to file
        with open('output/results.json', 'w') as f:
            json.dump(results, f, indent=2)
            
    except Exception as e:
        simulation_status['message'] = f'Error: {str(e)}'
        simulation_status['running'] = False

@app.route('/run-model', methods=['POST'])
def run_model():
    """Endpoint to start the GR4J simulation"""
    global simulation_status
    
    if simulation_status['running']:
        return jsonify({'error': 'Simulation already running'}), 400
    
    params = request.json
    print(f"Starting simulation with params: {params}")
    
    # Start simulation in background thread
    thread = threading.Thread(target=run_gr4j_simulation, args=(params,))
    thread.daemon = True
    thread.start()
    
    return jsonify({'message': 'Simulation started', 'id': 'current'})

@app.route('/status')
def get_status():
    """Endpoint to get simulation status"""
    return jsonify(simulation_status)

@app.route('/results')
def get_results():
    """Endpoint to get simulation results"""
    if simulation_status['results'] is None:
        return jsonify({'error': 'No results available'}), 404
    
    return jsonify(simulation_status['results'])

@app.route('/download/<file_type>')
def download_file(file_type):
    """Endpoint to download results files"""
    if file_type == 'csv':
        # Generate CSV from results
        results = simulation_status['results']
        if results:
            df = pd.DataFrame({
                'date': results['dates'],
                'observed_flow': results['observed_flow'],
                'simulated_flow': results['simulated_flow'],
                'simulated_volume': results['simulated_volume']
            })
            csv_path = 'output/results.csv'
            df.to_csv(csv_path, index=False)
            return send_file(csv_path, as_attachment=True)
    
    return jsonify({'error': 'File not found'}), 404

if __name__ == '__main__':
    os.makedirs('output', exist_ok=True)
    app.run(host='0.0.0.0', port=5000, debug=False)
