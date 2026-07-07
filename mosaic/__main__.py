"""Enable `python -m mosaic ...` as an alias for the `mosaic` CLI."""
from .cli import cli

if __name__ == "__main__":
    cli()
