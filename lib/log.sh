# lib/log.sh — bash logging helpers for the metabo.figures pipeline.
#
# Sourced by every bash entrypoint. Writes to $METABO_FIGURES_LOG using
# the unified line format documented in
# specs/001-figures-pipeline/contracts/log-format.md:
#
#   <ISO-8601-with-tz> [bash][<level>][<component>] <message>
#   <ISO-8601-with-tz> [LaTeX][<level>][<component>] <captured xelatex line>
#
# All writes go through printf with single line atomicity (POSIX
# O_APPEND, lines well under PIPE_BUF) so concurrent --jobs > 1
# invocations interleave safely.

: "${METABO_FIGURES_LOG:=}"

if [ -z "$METABO_FIGURES_LOG" ]; then
    METABO_FIGURES_LOG="logs/orphan-$$.log"
    mkdir -p "$(dirname "$METABO_FIGURES_LOG")"
    printf '[lib/log.sh] WARNING: METABO_FIGURES_LOG not set; using %s\n' \
        "$METABO_FIGURES_LOG" >&2
fi
export METABO_FIGURES_LOG

# log_line LEVEL COMPONENT MESSAGE...
#
# Emit a single bash-source line into the pipeline log.
log_line() {
    local level=$1
    local component=$2
    shift 2
    printf '%s [bash][%s][%s] %s\n' \
        "$(date -Iseconds)" "$level" "$component" "$*" \
        >> "$METABO_FIGURES_LOG"
}

# log_xelatex_capture COMPONENT CMD [ARGS...]
#
# Run a xelatex (or any latex-tooling) command, capture its stdout +
# stderr line-by-line, and prepend the unified-format prefix as
# [LaTeX][INFO][<component>] before appending to the pipeline log. Also
# preserves the command's exit status. Lines longer than 4 KiB are
# truncated with a trailing ellipsis (POSIX atomicity envelope).
log_xelatex_capture() {
    local component=$1
    shift
    local rc

    "$@" 2>&1 | while IFS= read -r line; do
        if [ "${#line}" -gt 4090 ]; then
            line="${line:0:4090}…"
        fi
        printf '%s [LaTeX][INFO][%s] %s\n' \
            "$(date -Iseconds)" "$component" "$line" \
            >> "$METABO_FIGURES_LOG"
    done
    rc=${PIPESTATUS[0]}

    if [ "$rc" -ne 0 ]; then
        log_line ERROR "$component" "xelatex exited with status $rc"
    fi
    return "$rc"
}
