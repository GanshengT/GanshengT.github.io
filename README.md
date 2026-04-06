Gansheng TAN personal web page was forked (then detached) by [Stuart Geiger](https://github.com/staeiou) from the [Minimal Mistakes Jekyll Theme](https://mmistakes.github.io/minimal-mistakes/), which is © 2016 Michael Rose and released under the MIT License. See LICENSE.md.

## Local Preview Cheatsheet

Run everything from:

```bash
cd "/Users/ganshengtan/Library/CloudStorage/Box-Box/Washu/projects/github_project/academicpages.github.io"
```

Start a local Jekyll preview:

```bash
bundle exec jekyll serve --livereload
```

Open:

```text
http://127.0.0.1:4000
```

Stop the server:

```bash
Ctrl-C
```

If the build gets weird:

```bash
bundle exec jekyll clean
bundle exec jekyll serve --livereload
```

If you change `_config.yml`, stop and restart the server.
