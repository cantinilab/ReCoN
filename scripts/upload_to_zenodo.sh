#!/bin/bash
#
# Upload ReCoN tutorial data to Zenodo
#
# Usage:
#   1. Export a token from https://zenodo.org/account/settings/applications/
#   2. Create a new deposit: ./scripts/upload_to_zenodo.sh
#   3. Extend an existing record: ./scripts/upload_to_zenodo.sh --new-version RECORD_ID
#
# For testing with sandbox:
#   ./scripts/upload_to_zenodo.sh --sandbox

set -e

# Configuration
DATA_DIR="data"
FILES=(
    "perturbation_tuto/rna.h5ad"
    "perturbation_tuto/rna_treated.h5ad"
    "perturbation_tuto/grn.csv"
    "build_grn_tuto/pbmc10x.h5mu"
)

SPATIAL_FILES=(
    "spatial_tuto/reference_heart_4k.h5ad"
    "spatial_tuto/visium18_all_spots_4k_genes.h5ad"
    "spatial_tuto/heart_grn_top250k.csv"
)

# Zenodo metadata
TITLE="ReCoN Tutorial Data"
DESCRIPTION="Tutorial data for the ReCoN (Reconstruction of Multicellular Systems) Python package. Includes single-cell, spatial transcriptomics, cell2location, and gene regulatory network inputs."
CREATOR_NAME="Trimbour, Rémi"
CREATOR_AFFILIATION="Institut Pasteur"
LICENSE="cc-by-4.0"
UPLOAD_TYPE="dataset"

# Parse arguments. A positional token remains supported for compatibility, but
# ZENODO_TOKEN is preferred because it does not put the token in shell history.
TOKEN="${ZENODO_TOKEN:-}"
SANDBOX=false
NEW_VERSION_ID=""
while [ "$#" -gt 0 ]; do
    case "$1" in
        --sandbox)
            SANDBOX=true
            shift
            ;;
        --new-version)
            NEW_VERSION_ID="${2:-}"
            if [ -z "$NEW_VERSION_ID" ]; then
                echo "ERROR: --new-version requires the latest record ID"
                exit 1
            fi
            shift 2
            ;;
        *)
            if [ -z "$TOKEN" ]; then
                TOKEN="$1"
                shift
            else
                echo "ERROR: Unknown argument: $1"
                exit 1
            fi
            ;;
    esac
done

if [ -z "$TOKEN" ]; then
    echo "ERROR: Set ZENODO_TOKEN or pass a token as the first argument."
    exit 1
fi

if [ "$SANDBOX" = true ]; then
    ZENODO_URL="https://sandbox.zenodo.org/api"
    echo "Using Zenodo SANDBOX"
else
    ZENODO_URL="https://zenodo.org/api"
    echo "Using Zenodo PRODUCTION"
fi

if [ -n "$NEW_VERSION_ID" ]; then
    UPLOAD_FILES=("${SPATIAL_FILES[@]}")
else
    UPLOAD_FILES=("${FILES[@]}" "${SPATIAL_FILES[@]}")
fi

# Check files exist
echo ""
echo "Checking files..."
for file in "${UPLOAD_FILES[@]}"; do
    filepath="$DATA_DIR/$file"
    if [ ! -f "$filepath" ]; then
        echo "  ERROR: $filepath not found"
        exit 1
    fi
    size=$(du -h "$filepath" | cut -f1)
    echo "  ✓ $file ($size)"
done

# Create a new deposit, or clone the latest published record into a new draft.
echo ""
if [ -n "$NEW_VERSION_ID" ]; then
    echo "1. Creating a new version of record $NEW_VERSION_ID..."
    RESPONSE=$(curl -fsS -X POST \
        "$ZENODO_URL/deposit/depositions/$NEW_VERSION_ID/actions/newversion" \
        -H "Authorization: Bearer $TOKEN")
    DRAFT_URL=$(echo "$RESPONSE" | jq -r '.links.latest_draft // empty')
    if [ -z "$DRAFT_URL" ]; then
        echo "ERROR: Zenodo did not return links.latest_draft"
        echo "$RESPONSE" | jq '{message, status, errors}'
        exit 1
    fi
    RESPONSE=$(curl -fsS "$DRAFT_URL" -H "Authorization: Bearer $TOKEN")
else
    echo "1. Creating new deposit..."
    RESPONSE=$(curl -fsS -X POST "$ZENODO_URL/deposit/depositions" \
        -H "Authorization: Bearer $TOKEN" \
        -H "Content-Type: application/json" \
        -d '{}')
fi

DEPOSIT_ID=$(echo "$RESPONSE" | jq -r '.id // empty')
BUCKET_URL=$(echo "$RESPONSE" | jq -r '.links.bucket // empty')

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

for file in "${UPLOAD_FILES[@]}"; do
    filepath="$DATA_DIR/$file"
    # Use basename for Zenodo URL (it doesn't support subdirectories)
    basename=$(basename "$file")
    size=$(du -h "$filepath" | cut -f1)
    echo ""
    echo "  Uploading $file ($size)..."
    
    UPLOAD_RESPONSE=$(curl -fsS -X PUT "$BUCKET_URL/$basename" \
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

# A new-version draft already inherits the complete metadata of the published
# record. Do not replace it: doing so would remove contributors and other
# fields that are not represented in the minimal metadata template below.
echo ""
if [ -n "$NEW_VERSION_ID" ]; then
    echo "3. Preserving inherited metadata for the new-version draft."
else
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
fi

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

if [ "$SANDBOX" = true ]; then
    DOWNLOAD_URL="https://sandbox.zenodo.org/records/$DEPOSIT_ID/files/"
else
    DOWNLOAD_URL="https://zenodo.org/records/$DEPOSIT_ID/files/"
fi

echo "TUTORIAL_DATA_URL = \"$DOWNLOAD_URL\""
echo "TUTORIAL_DATA_REGISTRY = {"
for file in "${UPLOAD_FILES[@]}"; do
    echo "    \"$file\": \"${FILE_HASHES[$file]}\","
done
echo "}"
echo ""
