# AMR-Intel Cloud Deployment

This is a Docker/Linux redesign of AMR-Intel. It removes Windows/WSL-specific commands and runs bioinformatics tools directly inside the cloud container.

## Files

- `app_cloud.py`: Cloud-compatible Streamlit app.
- `environment.yml`: Conda/Mamba bioinformatics environment.
- `Dockerfile`: Docker build file for full deployment.

## Recommended platforms

Use a Docker-supporting platform, such as Hugging Face Spaces Docker, Render Docker Web Service, Railway, Google Cloud Run, AWS ECS/EC2, or Azure Container Apps.

Streamlit Community Cloud is not recommended for the full genome-analysis version because it does not support custom Dockerfiles and the backend is heavy.

## Local Docker test

```bash
docker build -t amr-intel-cloud .
docker run -p 8501:8501 amr-intel-cloud
```

Open:

```text
http://localhost:8501
```

## Important

The Dockerfile preloads CARD, AMRFinderPlus, and ABRicate databases during build. The build can take time.
