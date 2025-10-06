import matplotlib.pyplot as plt
import numpy as np
from matplotlib.patches import FancyArrowPatch, Rectangle
import matplotlib.patheffects as path_effects

# Set up the figure with high resolution for quality
plt.figure(figsize=(14, 10), dpi=300)
ax = plt.gca()
ax.set_xlim(0, 100)
ax.set_ylim(0, 70)
ax.axis('off')

# Define gradient backgrounds for boxes
def gradient_fill(x, y, fill_color, ax=None, alpha=0.5):
    if ax is None:
        ax = plt.gca()
    
    # Create gradient
    z = np.empty((100, 100))
    for i in range(100):
        z[i, :] = i
    
    # Plot gradient
    rect = Rectangle((x, y), 20, 8, linewidth=2, edgecolor='black', facecolor='none')
    ax.add_patch(rect)
    
    # Add gradient within rectangle bounds
    gradient = ax.imshow(z, cmap=plt.cm.get_cmap('Blues'), 
                        extent=[x, x+20, y, y+8], 
                        aspect='auto', alpha=alpha, origin='lower')
    gradient.set_clip_path(rect)
    
    return rect

# Define vibrant colors
colors = {
    'bg': '#f0f8ff',  # Light blue background
    'arrow': '#FF6B6B',  # Coral red
    'box1': '#4ECDC4',  # Teal
    'box2': '#FF9F1C',  # Orange
    'box3': '#7A28CB',  # Purple
    'box4': '#2A9D8F',  # Green
    'box5': '#E76F51',  # Terracotta
    'box6': '#457B9D',  # Steel blue
    'title': '#1D3557'  # Dark blue
}

# Set the background color
fig = plt.gcf()
fig.patch.set_facecolor(colors['bg'])

# Create a beautiful title with shadow effect
title = plt.text(50, 65, 'EGFR Allosteric Binder Discovery Pipeline', 
                 fontsize=28, ha='center', fontweight='bold', color=colors['title'])
title.set_path_effects([path_effects.withStroke(linewidth=3, foreground='white')])

subtitle = plt.text(50, 60, 'Computational workflow for identifying and evaluating allosteric binding sites', 
                    fontsize=14, ha='center', color='#333333')

# Function to draw fancy boxes with shadows and gradients
def draw_step(x, y, title, description, color, step_number):
    # Add shadow
    shadow = Rectangle((x+0.5, y-0.5), 20, 8, linewidth=0, 
                       facecolor='gray', alpha=0.3)
    ax.add_patch(shadow)
    
    # Add box with gradient
    rect = gradient_fill(x, y, color, ax, alpha=0.7)
    
    # Add step circle
    circle = plt.Circle((x-3, y+4), 2.5, color=color, alpha=0.9, 
                      edgecolor='black', linewidth=1.5)
    ax.add_patch(circle)
    
    # Add step number
    step_text = plt.text(x-3, y+4, str(step_number), ha='center', va='center', 
                         fontsize=14, fontweight='bold', color='white')
    
    # Add titles and descriptions
    plt.text(x+10, y+6, title, ha='center', fontsize=12, 
             fontweight='bold', color='#1D3557')
    plt.text(x+10, y+3, description, ha='center', fontsize=9, 
             color='#333333', wrap=True)

# Draw workflow steps in a zigzag pattern
steps = [
    {
        'title': 'Input Preparation',
        'desc': 'Process MMCIF files\nfor structural analysis',
        'color': colors['box1']
    },
    {
        'title': 'Sample Building',
        'desc': 'Compile IC50 data\nand sample information',
        'color': colors['box2']
    },
    {
        'title': 'Site Detection',
        'desc': 'Unsupervised detection\nof allosteric sites',
        'color': colors['box3']
    },
    {
        'title': 'Pose Selection',
        'desc': 'Unbiased selection of\noptimal binding poses',
        'color': colors['box4']
    },
    {
        'title': 'Contact Analysis',
        'desc': 'PLIP analysis of\nprotein-ligand interactions',
        'color': colors['box5']
    },
    {
        'title': 'Scoring & Clustering',
        'desc': 'Evaluate and group\npotential binders',
        'color': colors['box6']
    }
]

# Position the boxes in a zigzag flow
positions = [
    (10, 45),  # Step 1
    (40, 45),  # Step 2
    (70, 45),  # Step 3
    (70, 30),  # Step 4
    (40, 30),  # Step 5
    (10, 30),  # Step 6
]

# Draw all steps
for i, (step, pos) in enumerate(zip(steps, positions)):
    draw_step(pos[0], pos[1], step['title'], step['desc'], step['color'], i+1)

# Draw arrows connecting the steps
arrows = [
    ((30, 49), (40, 49)),    # Step 1 to 2
    ((60, 49), (70, 49)),    # Step 2 to 3
    ((70, 45), (70, 38)),    # Step 3 to 4
    ((70, 34), (60, 34)),    # Step 4 to 5
    ((40, 34), (30, 34)),    # Step 5 to 6
]

for start, end in arrows:
    arrow = FancyArrowPatch(start, end, 
                          connectionstyle='arc3,rad=0.1',
                          arrowstyle='fancy,head_width=4,head_length=6',
                          facecolor=colors['arrow'], linewidth=2,
                          edgecolor='#333333', alpha=0.8)
    ax.add_patch(arrow)

# Add a data flow section at the bottom
plt.text(50, 20, 'Data Flow', ha='center', fontsize=16, 
        fontweight='bold', color=colors['title'])

# Data structures
data_items = [
    ('MMCIF Files', 15, 12, '#82AAFF'),
    ('Standardized Data', 35, 12, '#C3E88D'),
    ('Binding Sites', 55, 12, '#F78C6C'),
    ('Ranked Molecules', 75, 12, '#BB80B3')
]

for name, x, y, color in data_items:
    rect = Rectangle((x-7, y-2), 14, 4, linewidth=1.5, 
                   edgecolor='black', facecolor=color, alpha=0.7)
    ax.add_patch(rect)
    plt.text(x, y, name, ha='center', va='center', fontsize=10, 
            fontweight='bold', color='#333333')

# Connect data items with dotted arrows
for i in range(len(data_items)-1):
    arrow = FancyArrowPatch((data_items[i][1]+7, data_items[i][2]), 
                          (data_items[i+1][1]-7, data_items[i+1][2]),
                          connectionstyle='arc3,rad=0',
                          arrowstyle='->,head_width=1.5,head_length=2',
                          linestyle='dotted', linewidth=1.5,
                          color='#333333')
    ax.add_patch(arrow)

# Add a credit line
plt.text(50, 5, 'Created for boltz2-egfr-allosteric-binder-discovery project',
        fontsize=9, color='gray', ha='center')

# Save the image in high resolution
plt.savefig('workflow_diagram.png', dpi=300, bbox_inches='tight', 
           facecolor=fig.get_facecolor())
plt.savefig('workflow_diagram.pdf', format='pdf', bbox_inches='tight',
           facecolor=fig.get_facecolor())

print("Workflow diagram created successfully!")
print("Files saved as 'workflow_diagram.png' and 'workflow_diagram.pdf'")
