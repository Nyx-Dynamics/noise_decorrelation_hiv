from pathlib import Path

# Get project root
base_path = Path(__file__).resolve().parent.parent
data_path = base_path / 'data'

print(f"Project root: {base_path}")
print(f"Data directory: {data_path}")
print(f"Data directory exists: {data_path.exists()}")

if data_path.exists():
    print("\nContents of data/:")
    for item in sorted(data_path.iterdir()):
        print(f"  {item.name}")

    # Check subdirectories
    for subdir in ['individual', 'extracted', 'raw']:
        sub_path = data_path / subdir
        if sub_path.exists():
            print(f"\nContents of data/{subdir}/:")
            for item in sorted(sub_path.iterdir()):
                print(f"  {item.name}")