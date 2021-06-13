FROM julia:1.6
WORKDIR app
COPY . .
RUN julia --project -e 'using Pkg; Pkg.instantiate()'
ENTRYPOINT julia --project -e 'using DoubleTrouble; main()'
