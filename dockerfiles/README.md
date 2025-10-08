# Dockerfiles

Custom Docker image definitions for tools used in this project.

Container definitions are written using the Dockerfile format to allow for portability
and cross-platform compatibity ( i.e. that the container will run on the majority of
container platforms such as Docker, Singularity/Apptainer, Podman, Shifter, CharlieCloud, etc).

## Building a container

Containers are built locally and then pushed to `ghcr.io` (the GitHub Container Registry) where they're hosted
and linked to this repository. Docker needs to be installed locally to build and push. 

> **Note**
> To push to the GitHub Container Registry, you need to first provide Docker with credentials from GitHub.
>
> Use a GitHub Personal Access Token (PAT) to authenticate yourself. In your web-browser, click on your
user-icon in GitHub in the top-right of the screen. Then go to `Settings > Developer Settings > Personal Access Tokens`,
click on “Generate new token”, and enable “write:packages”. Provide a note, duration, and enable other settings
as necessary, followed by selecting “Generate token”. Copy the token.
>
> Then in your terminal:
>
> ```bash
> GHCR_TOKEN=<your_token>
> echo "$GHCR_TOKEN" | docker login ghcr.io -u <GitHub username> --password-stdin
> ```

Check the Dockerfile for any build options (`ARG`) which need to be provided. Then
build the image in the root of this project with:

```bash
docker build --platform=linux/amd64 -t ghcr.io/nbisweden/<image>:<tag> -f </path/to/Dockerfile> .
```

where:
- `--platform=linux/amd64` means the image will be built to run on Intel chips (amd64).
- `<image>` is the name of the software.
- `<tag>` is the either the software version number if building from a versioned release, or a commit id if building from a version
controlled source without release versions.
- `</path/to/Dockerfile>` is the path to the Dockerfile, e.g. `./dockerfiles/fastk.Dockerfile`.

Once the image has been built (and you've authenticated to ghcr via docker), you can push the image to `ghcr.io`.

```bash
docker push ghcr.io/nbisweden/<image>:<tag>
```