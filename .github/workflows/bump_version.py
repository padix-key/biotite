import argparse
import json
from pathlib import Path
import re


BIOTITE_URL = "https://www.biotite-python.org"
SEMVER_REGEX = r"\d+\.\d+\.\d+"


def is_valid_semantic_version(version):
    return re.match(f"^{SEMVER_REGEX}$", version) is not None


def update_package_version(version, package_root):
    init_file_path = package_root / "__init__.py"
    with open(init_file_path, "r") as init_file:
        content = init_file.read()

    updated_content = re.sub(SEMVER_REGEX, version, content)

    with open(init_file_path, "w") as init_file:
        init_file.write(updated_content)


def add_version_to_switcher(version, switcher_file):
    with open(switcher_file, "r") as file:
        version_config = json.load(file)

    version_wo_patch = ".".join(version.split(".")[:2])
    # Remove existing entry for this minor.major version, if existing
    version_config = [
        entry for entry in version_config if entry["name"] != version_wo_patch
    ]
    # Remove 'preferred' tag from existing entries,
    # as the new version is the new preferred version
    for entry in version_config:
        if "preferred" in entry:
            del entry["preferred"]
    # Add entry for the new version
    version_config.append({
        "name": version_wo_patch,
        "version": version,
        "url": f"{BIOTITE_URL}/{version}/",
        "preferred": True
    })

    with open(switcher_file, "w") as file:
        json.dump(version_config, file, indent=4)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Create a new version by updating the necessary files",
    )
    parser.add_argument("version", type=str, help="The version to be created")
    args = parser.parse_args()
    version = args.version
    if not is_valid_semantic_version(version):
        raise ValueError(f"{version} is not a valid version string")

    project_root = Path(__file__).parents[2]
    update_package_version(
        version,
        project_root / "src" / "biotite"
    )
    add_version_to_switcher(
        version,
        project_root / "doc" / "static" / "switcher.json"
    )