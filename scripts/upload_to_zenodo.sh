#!/bin/bash
#
# Upload ReCoN tutorial data to Zenodo
#
# Usage:
#   1. Get your Zenodo API token from: https://zenodo.org/account/settings/applications/
#   2. Run: ./scripts/upload_to_zenodo.sh <ZENODO_TOKEN>
#
# For testing with sandbox:
#   ./scripts/upload_to_zenodo.sh <SANDBOX_TOKEN> --sandbox

set -e

# Configuration
DATA_DIR="data"
FILES=(
    "perturbation_tuto/rna.h5ad"
    "perturbation_tuto/rna_treated.h5ad"
    "perturbation_tuto/grn.csv"
    "build_grn_tuto/pbmc10x.h5mu"
)

# Zenodo metadata
TITLE="ReCoN Tutorial Data - Single-cell RNA-seq and GRN"
DESCRIPTION="Tutorial data for the ReCoN (Reconstruction of Multicellular Systems) Python package. Includes scRNA-seq data and pre-computed gene regulatory network."
CREATOR_NAME="Trimbour, Rémi"
CREATOR_AFFILIATION="Institut Pasteur"
LICENSE="cc-by-4.0"
UPLOAD_TYPE="dataset"

# Parse arguments
if [ -z "$1" ]; then
    echo "Usage: $0 <ZENODO_TOKEN> [--sandbox]"
    echo ""
    echo "Get your token from: https://zenodo.org/account/settings/applications/"
    exit 1
fi

TOKEN="$1"

if [ "$2" == "--sandbox" ]; then
    ZENODO_URL="https://sandbox.zenodo.org/api"
    echo "Using Zenodo SANDBOX"
else
    ZENODO_URL="https://zenodo.org/api"
    echo "Using Zenodo PRODUCTION"
fi

# Check files exist
echo ""
echo "Checking files..."
for file in "${FILES[@]}"; do
    filepath="$DATA_DIR/$file"
    if [ ! -f "$filepath" ]; then
        echo "  ERROR: $filepath not found"
        exit 1
    fi
    size=$(du -h "$filepath" | cut -f1)
    echo "  ✓ $file ($size)"
done

# Create new deposit
echo ""
echo "1. Creating new deposit..."
RESPONSE=$(curl -s -X POST "$ZENODO_URL/deposit/depositions" \
    -H "Authorization: Bearer $TOKEN" \
    -H "Content-Type: application/json" \
    -d '{}')

DEPOSIT_ID=$(echo "$RESPONSE" | grep -o '"id": [0-9]*' | head -1 | grep -o '[0-9]*')
BUCKET_URL=$(echo "$RESPONSE" | grep -o '"bucket": "[^"]*"' | head -1 | cut -d'"' -f4)

if [ -z "$DEPOSIT_ID" ]; then
    echo "  ERROR: Failed to create deposit"
    echo "$RESPONSE"
    exit 1
fi

echo "  Deposit ID: $DEPOSIT_ID"
echo "  Bucket URL: $BUCKET_URL"

# Upload files
echo ""
echo "2. Uploading files..."
declare -A FILE_HASHES

for file in "${FILES[@]}"; do
    filepath="$DATA_DIR/$file"
    # Use basename for Zenodo URL (it doesn't support subdirectories)
    basename=$(basename "$file")
    size=$(du -h "$filepath" | cut -f1)
    echo ""
    echo "  Uploading $file ($size)..."
    
    UPLOAD_RESPONSE=$(curl -s -X PUT "$BUCKET_URL/$basename" \
        -H "Authorization: Bearer $TOKEN" \
        -H "Content-Type: application/octet-stream" \
        --data-binary @"$filepath")
    
    # Check for errors
    if echo "$UPLOAD_RESPONSE" | grep -q '"status":\s*4[0-9][0-9]\|"status":\s*5[0-9][0-9]'; then
        echo "  ERROR uploading $file:"
        echo "$UPLOAD_RESPONSE"
        exit 1
    fi
    
    # Calculate SHA256
    hash=$(sha256sum "$filepath" | cut -d' ' -f1)
    FILE_HASHES[$file]="sha256:$hash"
    echo "  ✓ Uploaded. SHA256: sha256:$hash"
done

# Update metadata
echo ""
echo "3. Updating metadata..."
METADATA=$(cat <<EOF
{
    "metadata": {
        "title": "$TITLE",
        "upload_type": "$UPLOAD_TYPE",
        "description": "$DESCRIPTION",
        "creators": [{"name": "$CREATOR_NAME", "affiliation": "$CREATOR_AFFILIATION"}],
        "license": "$LICENSE",
        "related_identifiers": [
            {"identifier": "https://github.com/cantinilab/ReCoN", "relation": "isSupplementTo", "scheme": "url"}
        ]
    }
}
EOF
)

curl -s -X PUT "$ZENODO_URL/deposit/depositions/$DEPOSIT_ID" \
    -H "Authorization: Bearer $TOKEN" \
    -H "Content-Type: application/json" \
    -d "$METADATA" > /dev/null

echo "  ✓ Metadata updated"

# Print summary
echo ""
echo "========================================"
echo "DEPOSIT READY"
echo "========================================"
echo ""
echo "Deposit ID: $DEPOSIT_ID"
echo "View/edit at: ${ZENODO_URL/api/}/deposit/$DEPOSIT_ID"
echo ""
echo "To publish, run:"
echo "  curl -X POST '$ZENODO_URL/deposit/depositions/$DEPOSIT_ID/actions/publish' -H 'Authorization: Bearer $TOKEN'"
echo ""
echo "========================================"
echo "UPDATE load_data.py WITH:"
echo "========================================"
echo ""

if [ "$2" == "--sandbox" ]; then
    DOWNLOAD_URL="https://sandbox.zenodo.org/records/$DEPOSIT_ID/files/"
else
    DOWNLOAD_URL="https://zenodo.org/records/$DEPOSIT_ID/files/"
fi

echo "TUTORIAL_DATA_URL = \"$DOWNLOAD_URL\""
echo "TUTORIAL_DATA_REGISTRY = {"
for file in "${FILES[@]}"; do
    echo "    \"$file\": \"${FILE_HASHES[$file]}\","
done
echo "}"
echo ""
