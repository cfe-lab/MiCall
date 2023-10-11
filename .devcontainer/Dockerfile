FROM cfelab/micall:v7.15.13

# Uninstall some Python packages that don't upgrade correctly.
RUN python -m pip install --upgrade pip && \
    python -m pip uninstall -y python-Levenshtein Levenshtein genetracks

ENTRYPOINT bash