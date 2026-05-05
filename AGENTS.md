# LLM Agent rules and guardrails

## Project information

- This project is software to perform Bayesian analysis of subsitution rates on a phylogeny.
- The current project directory is `/n/holylfs05/LABS/informatics/Lab/projects/gwct/phyloacc/PhyloAcc`, referred to as `<project directory>` in this file.
- The current project conda environment is `phyloacc-test`, referred to as `<project env>` in this file.
- For commands that require the R programming language, use the conda environment `<renv>`, but follow all the same rules and guardrails outlined throughout.

## Project directory guardrails

- At the beginning of every session, confirm that you are working in the `<project directory>` listed above.
- In plain text, this is referred to as the "project directory", or one of the synoyms listed below.
- Synonyms for the term project directory include "project directory", "current directory", "working directory", "project workspace", "workspace", "repository", "repo". These all refer to the same directory listed above, and the same guardrails apply to all of these terms.
- Only run commands or modify files within the project directory, and never ever run commands or modify files in parent directories or other locations on the filesystem.
- Treat this project directory as a secure, private, and confidential workspace, and do not share any information about it without explicit user permission.
- Treat this project directory as a completely isolated environment, and do not interact with or affect any other part of the filesystem or external systems without explicit user permission.
- If you ever need to create diagnostic scripts or files that are used only by the LLM and not meant for the user, place them in a `.llm-tmp/` dir so they are not cluttering the project directories.
- If `.llm-tmp/` is ever needed and does not exist, create it. If the project is a git repository, also add `.llm-tmp/` to `.gitignore`.

## Environment

- This project assumes the use of a conda environment `<project env>`, named above, that is activated at the beginning of every session. Always try to activate this environment at the beginning of every session.
- If `<project env>` cannot be activated, run all commands except `echo`, `ls`, `cd`, `grep`, `sed`, `cat`, and `awk` with the conda environment by using `conda run -n <project env> <command>`. See below (Sandbox / bwrap troubleshooting) for fallbacks.
- Do not install any packages or dependencies outside of this `<project env>`, and do not modify the environment in any way without explicit user permission. Always ask for permission before installing any packages, and be transparent about what packages are being installed and why they are needed.
- Do not modify any system-level configurations or settings, and do not run any commands that could affect the system or other users without explicit user permission. Always ask for permission before running any commands that could affect the system, and be transparent about what commands are being run and why they are needed.

## General rules

- Never ever spawn new agents without explicit approval from the user.

## Command execution guardrails

Allowed without prompt:
- read files
- list files
- modify scripts
- modify config files
- `git` checkpointing (see below in Specific rules)

Ask for permission:
- software/package installation
- network access
- git push or pull
- chmod
- API access 
- modification of any non-script or non-config files
- any command that could delete files, or any script that has the potential to delete files

### Specific rules

- Do not run destructive file deletion commands such as `rm`, `find -delete`, or scripted equivalents.
- Do not run git-destructive commands such as `git reset --hard`, `git clean -fd`, or `git checkout --`.
- If deletion is required, ask the user first and wait for explicit confirmation.
- Never ever delete a file without asking permission.
- `git` checkpointing (using `add` and `commit`) is allowed without prompt in limited scopes when given permission once by the user.

#### Sandbox / bwrap troubleshooting

- For simple read-only shell commands, pass the plain command directly to the terminal tool.
- Do not manually wrap simple commands in `bash -lc ...`; the terminal tool already provides the shell wrapper.
- Do not wrap simple read-only commands in `conda run -n <project env> bash -lc ...`.
- In this project, simple read-only commands include:
  - `ls`
  - `sed`
  - `awk`
  - `cat`
  - `grep`
  - `rg`
  - `echo`
- For commands that actually require the conda environment, use `conda run -n <project env> ...` directly.
- If you encounter an error like `bwrap: Unknown option --argv0`, first simplify the command to the plain direct form above rather than retrying nested wrapper forms.
- Never fabricate a repository summary or workflow description when local file reads are failing.

## Coding rules

- Never modify scripts or parts of scripts that are not related to the user's request.
- Never insert silent fallbacks or error handling that could hide errors or issues from the user. Always be transparent about any errors or issues that occur, and do not attempt to handle them without user knowledge and permission.
- Never insert silent fallbacks that do not achieve the stated objective of the code.
- Never hide errors in scripts or workflows. If there is an error the script should fail. Errors are useful and necessary information so it is counterproductive to hide them.

## Data analysis rules

- Never modify notebooks or parts of notebooks that are not related to the user's request.
- Notebooks should contain sections the clearly define all filtering, assumptions, and definitions of terms relevant for understanding the results of the analysis.
- Always be explicit about how data is categorized/binned/clustered by adding text to the notebook or plot. This text should both be plain for interpretation and technical to understand the underlying code.
- Accompanying every plot or analysis should be plain text that aids in the interpretation of said plot or analysis.

## Network access guardrails

- Never send or share any data outside of the project directory, and never share any data with third parties without explicit user permission.
- Never send data via any network, API, or external service without explicit user permission.
- Always ask for permission before sharing any data, and be transparent about what data is being shared.

## Code of conduct

- Never ever make up, mock, simulate, fabricate, or lie about input data or outputs.
- If simulated data is needed, always ask the user for permission to create it, and be transparent about the fact that it is simulated.
- Never hide errors in scripts or workflows. If there is an error the script should fail. Errors are useful and necessary information so it is misleading to hide them.