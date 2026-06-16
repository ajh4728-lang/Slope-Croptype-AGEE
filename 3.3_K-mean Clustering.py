"""
K-means clustering of natural-rainfall runoff events based on TN and TP loads.

This script reproduces the cluster analysis reported in the manuscript
"Crop Type and Slope Interactions Govern Nutrient Loss Patterns in South
Korean Upland Fields: An Integrated Field Monitoring Analysis."

Workflow:
    1. Load the compiled monitoring dataset.
    2. Filter to Rainfall <= 200 mm/event and SlopeGradient <= 30%.
    3. Standardize TN and TP loads.
    4. Evaluate cluster validity (Elbow / distance-to-line and silhouette).
    5. Fit the final K-means model and relabel clusters by ascending mean TP.
    6. Export summary tables and figures.

Requirements:
    pandas, numpy, scikit-learn, matplotlib, scipy
"""

import pandas as pd
import numpy as np
from sklearn.cluster import KMeans
from sklearn.preprocessing import StandardScaler
from sklearn.metrics import silhouette_score
import matplotlib.pyplot as plt
from matplotlib import rcParams
from scipy.spatial import ConvexHull

# 1. Load the dataset
file_path = 'metadata.xlsx'
data = pd.read_excel(file_path)

# Set Arial as the default font
rcParams['font.family'] = 'Arial'

# 2. Keep only the required columns and drop rows with missing values
data_cleaned = data.dropna(
    subset=['CropType', 'SlopeGradient', 'TP', 'TN', 'Rainfall']
).copy()

# 3. Filter to Rainfall <= 200 mm/event and SlopeGradient <= 30%
data_filtered = data_cleaned[
    (data_cleaned['Rainfall'] <= 200) &
    (data_cleaned['SlopeGradient'] <= 30)
].copy()

# 4. Select clustering variables (TN, TP)
selected_data = data_filtered[['TN', 'TP']].copy()

# 5. Standardize the variables (zero mean, unit variance)
scaler = StandardScaler()
scaled_data = scaler.fit_transform(selected_data)

# 6. Compute WSS, distance-to-line, and silhouette coefficient across k
max_k = min(15, len(scaled_data) - 1)  # limit to n_samples - 1 for silhouette
k_range = range(1, max_k + 1)

wss_list = []
metric_rows = []

for k in k_range:
    kmeans = KMeans(
        n_clusters=k,
        init='k-means++',
        random_state=23, n_init=10
    )

    labels = kmeans.fit_predict(scaled_data)
    wss = kmeans.inertia_
    wss_list.append(wss)

    row = {
        'k': k,
        'WSS': wss,
        'Silhouette_score': np.nan
    }

    # The silhouette coefficient is defined only for k >= 2
    if k >= 2:
        row['Silhouette_score'] = silhouette_score(scaled_data, labels)

    metric_rows.append(row)

metrics_df = pd.DataFrame(metric_rows)


# 7. Automatic elbow-point detection
#    The elbow is the point on the normalized WSS curve with the maximum
#    perpendicular distance to the line connecting the first and last points.
def compute_elbow_distance_table(k_values, wss_values):
    k = np.array(k_values, dtype=float)
    w = np.array(wss_values, dtype=float)

    if len(k) < 3:
        elbow_k = int(k[np.argmin(w)])

        distance_df = pd.DataFrame({
            'k': k.astype(int),
            'WSS': w,
            'k_normalized': np.nan,
            'WSS_normalized': np.nan,
            'Distance_to_line': np.nan,
            'Distance_rank': np.nan,
            'Distance_percent_of_max': np.nan
        })

        return elbow_k, distance_df

    # Min-max normalization to the [0, 1] range
    k_norm = (k - k.min()) / (k.max() - k.min())
    w_norm = (w - w.min()) / (w.max() - w.min())

    # Line connecting the first and last points
    p1 = np.array([k_norm[0], w_norm[0]])
    p2 = np.array([k_norm[-1], w_norm[-1]])

    line_vec = p2 - p1
    line_unit = line_vec / np.linalg.norm(line_vec)

    distances = []

    for i in range(len(k_norm)):
        p = np.array([k_norm[i], w_norm[i]])
        v = p - p1
        proj = np.dot(v, line_unit) * line_unit
        distance = np.linalg.norm(v - proj)
        distances.append(distance)

    distance_df = pd.DataFrame({
        'k': k.astype(int),
        'WSS': w,
        'k_normalized': k_norm,
        'WSS_normalized': w_norm,
        'Distance_to_line': distances
    })

    distance_df['Distance_rank'] = distance_df['Distance_to_line'].rank(
        ascending=False,
        method='min'
    ).astype(int)

    distance_df['Distance_percent_of_max'] = (
        distance_df['Distance_to_line'] /
        distance_df['Distance_to_line'].max() * 100
    )

    elbow_k = int(distance_df.loc[
        distance_df['Distance_to_line'].idxmax(), 'k'
    ])

    return elbow_k, distance_df


# 7-1. Compute the elbow point and the distance table
optimal_k, elbow_distance_df = compute_elbow_distance_table(
    list(k_range),
    wss_list
)

# 7-2. Merge the distance information into metrics_df
metrics_df = metrics_df.merge(
    elbow_distance_df[
        ['k', 'Distance_to_line', 'Distance_rank', 'Distance_percent_of_max']
    ],
    on='k',
    how='left'
)

# 7-3. Optimal k based on the silhouette coefficient
valid_metrics = metrics_df.dropna(subset=['Silhouette_score']).copy()

sil_optimal_k = int(valid_metrics.loc[
    valid_metrics['Silhouette_score'].idxmax(), 'k'
])

print("\n================ Cluster validity results ================")
print(f"Elbow method optimal k       : {optimal_k}")
print(f"Silhouette score optimal k   : {sil_optimal_k} "
      f"(score = {valid_metrics.loc[valid_metrics['k'] == sil_optimal_k, 'Silhouette_score'].iloc[0]:.3f})")

if optimal_k == sil_optimal_k:
    print("-> Both indices suggest the same k")
else:
    print("-> The two indices suggest different k; review together with domain interpretation")

print("\nFull metric table:")
print(metrics_df.round(4))

print("\nElbow distance-to-line table:")
print(elbow_distance_df.round(4))

# Save the validity tables
metrics_df.to_csv(
    'KMeans_cluster_validity_indices.csv',
    index=False,
    encoding='utf-8-sig'
)

elbow_distance_df.to_csv(
    'KMeans_elbow_distance_to_line.csv',
    index=False,
    encoding='utf-8-sig'
)

print("\nSaved -> KMeans_cluster_validity_indices.csv")
print("Saved -> KMeans_elbow_distance_to_line.csv")


# 8. Elbow plot (WSS vs k)
plt.figure(figsize=(10, 7), dpi=500)
plt.plot(
    metrics_df['k'],
    metrics_df['WSS'],
    marker='o',
    label='WSS',
    color='black',
    markersize=5,
    linewidth=1
)
plt.axvline(
    x=optimal_k,
    color='red',
    linestyle='--',
    label=f'Elbow k = {optimal_k}'
)
plt.xlabel('Number of Clusters (k)', fontsize=20)
plt.ylabel('Within-Cluster Sum of Squares (WSS)', fontsize=20)
plt.xticks(fontsize=14)
plt.yticks(fontsize=14)
plt.legend(fontsize=16)
plt.grid(False)
plt.tight_layout()
plt.savefig('Elbow_WSS_KMeans_plot_filtered.png')
plt.show()


# 8-1. Distance-to-line plot
plt.figure(figsize=(10, 7), dpi=500)
plt.plot(
    elbow_distance_df['k'],
    elbow_distance_df['Distance_to_line'],
    marker='o',
    label='Distance to line',
    color='black',
    markersize=5,
    linewidth=1
)
plt.axvline(
    x=optimal_k,
    color='red',
    linestyle='--',
    label=f'Max distance k = {optimal_k}'
)
plt.xlabel('Number of Clusters (k)', fontsize=20)
plt.ylabel('Normalized Distance to Line', fontsize=20)
plt.xticks(fontsize=14)
plt.yticks(fontsize=14)
plt.legend(fontsize=16)
plt.grid(False)
plt.tight_layout()
plt.savefig('Elbow_Distance_to_Line_KMeans_plot_filtered.png')
plt.show()


# 8-2. Silhouette plot
plt.figure(figsize=(10, 7), dpi=500)
plt.plot(
    valid_metrics['k'],
    valid_metrics['Silhouette_score'],
    marker='o',
    label='Average Silhouette',
    color='black',
    markersize=5,
    linewidth=1
)
plt.axvline(
    x=sil_optimal_k,
    color='red',
    linestyle='--',
    label=f'Optimal k = {sil_optimal_k}'
)
plt.xlabel('Number of Clusters (k)', fontsize=20)
plt.ylabel('Average Silhouette Coefficient', fontsize=20)
plt.xticks(fontsize=14)
plt.yticks(fontsize=14)
plt.legend(fontsize=16)
plt.grid(False)
plt.tight_layout()
plt.savefig('Silhouette_KMeans_plot_filtered.png')
plt.show()


# 9. Fit the final K-means model with the selected k
kmeans_final = KMeans(
    n_clusters=optimal_k,
    init='k-means++',
    random_state=23
)

clusters_raw = kmeans_final.fit_predict(scaled_data)


# 9-1. Relabel clusters by ascending mean TP load
#      Cluster 1 = lowest mean TP
#      Higher cluster numbers correspond to higher mean TP
tmp = selected_data.copy()
tmp['raw'] = clusters_raw

order = tmp.groupby('raw')['TP'].mean().sort_values(ascending=True).index
relabel = {old: new + 1 for new, old in enumerate(order)}
clusters = np.array([relabel[c] for c in clusters_raw])

print("\nraw cluster -> final cluster relabel mapping:")
print(relabel)


# 10. Merge results and export
selected_data = selected_data.reset_index(drop=True)
data_filtered = data_filtered.reset_index(drop=True)

selected_data['Cluster'] = clusters

result = pd.concat([data_filtered, selected_data['Cluster']], axis=1)

result.to_csv(
    'KMeans_clustering_result_filtered.csv',
    index=False,
    encoding='utf-8-sig'
)

print("Clustering result saved -> 'KMeans_clustering_result_filtered.csv'")


# 11. Compute per-cluster mean TN, TP, and slope gradient
summary = result.groupby('Cluster')[['TP', 'TN', 'SlopeGradient']].mean().round(2)
summary['Sample Count'] = result.groupby('Cluster').size()

# 12. Print the summary
print("\nPer-cluster mean TN/TP and slope gradient:")
print(summary.sort_values(by=['TP', 'TN'], ascending=False))

summary.to_csv(
    'KMeans_cluster_summary.csv',
    encoding='utf-8-sig'
)

print("Cluster summary saved -> 'KMeans_cluster_summary.csv'")


# 13. Scatter plot of clusters (TN vs TP) with convex hulls
markers = ['o', 's', '^', 'D', 'v', '*', 'P', 'X', '<', '>']
unique_clusters = sorted(result['Cluster'].unique())

# Assign colors automatically based on the number of clusters
palette = plt.cm.tab10(np.linspace(0, 1, 10))
cluster_colors = {cid: palette[i % 10] for i, cid in enumerate(unique_clusters)}

plt.figure(figsize=(10, 7), dpi=500)

for idx, cluster_id in enumerate(unique_clusters):
    cluster_data = result[result['Cluster'] == cluster_id]
    marker = markers[idx % len(markers)]
    color = cluster_colors[cluster_id]

    # Scatter points
    plt.scatter(
        cluster_data['TN'],
        cluster_data['TP'],
        label=f'Cluster {cluster_id}',
        alpha=0.8,
        marker=marker,
        color=color,
        edgecolor='black',
        linewidth=0.5
    )

    # Convex hull
    if len(cluster_data) >= 3:
        points = cluster_data[['TN', 'TP']].values
        hull = ConvexHull(points)
        hull_points = points[hull.vertices]

        # Fill the hull interior
        plt.fill(
            hull_points[:, 0],
            hull_points[:, 1],
            color=color,
            alpha=0.2,
            zorder=0,
            linewidth=0
        )

        # Hull outline
        hull_points = np.append(hull_points, [hull_points[0]], axis=0)

        plt.plot(
            hull_points[:, 0],
            hull_points[:, 1],
            linestyle='-',
            linewidth=1.5,
            color=color
        )

# Axis settings
plt.xlabel('TN Load (kg/ha/event)', fontsize=20)
plt.ylabel('TP Load (kg/ha/event)', fontsize=20)
plt.xticks(fontsize=14)
plt.yticks(fontsize=14)
plt.legend(fontsize=16)
plt.grid(False)
plt.tight_layout()
plt.savefig('Cluster_Scatter_TN_TP_ColoredHull.png')
plt.show()