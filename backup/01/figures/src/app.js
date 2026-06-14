const d3 = window.d3;

const state = {
  data: null,
  resourceMetric: "total_size_bytes",
  entityFacet: "entity_type",
  relationFacet: "predicate",
  annotationDomain: "entity"
};

const palette = {
  teal: "#2a8c82",
  blue: "#4169a8",
  amber: "#c6861a",
  coral: "#c95d4f",
  green: "#5e8f3d",
  violet: "#7560a8",
  gray: "#8a9690"
};

const categorical = d3.scaleOrdinal([
  palette.teal,
  palette.blue,
  palette.amber,
  palette.coral,
  palette.green,
  palette.violet,
  "#64756f",
  "#9a7b42"
]);

const taxonNames = new Map([
  ["9606", "Human"],
  ["10090", "Mouse"],
  ["10116", "Rat"],
  ["9913", "Cattle"],
  ["7955", "Zebrafish"],
  ["8364", "African clawed frog"],
  ["9598", "Chimpanzee"],
  ["9615", "Dog"],
  ["9544", "Rhesus macaque"],
  ["9031", "Chicken"],
  ["4932", "Yeast"],
  ["3702", "Arabidopsis"],
  ["6239", "C. elegans"],
  ["7227", "Fruit fly"],
  ["562", "E. coli"]
]);

const metricLabels = {
  total_size_bytes: "Storage",
  entity_count: "Entities",
  association_count: "Associations",
  interaction_count: "Interactions",
  identifier_count: "Identifiers",
  ontology_term_count: "Ontology terms"
};

const landscapeRows = [
  {
    group: "Relations",
    name: "Associations",
    metric: "association_count",
    color: "#d49a2a"
  },
  {
    group: "Relations",
    name: "Interactions",
    metric: "interaction_count",
    color: "#12a9b6"
  },
  {
    group: "Entities",
    name: "Entity records",
    metric: "entity_count",
    color: "#4677b8"
  },
  {
    group: "Entities",
    name: "Identifiers",
    metric: "identifier_count",
    color: "#7a62b5"
  },
  {
    group: "Ontology",
    name: "Ontology terms",
    metric: "ontology_term_count",
    color: "#63a64e"
  },
  {
    group: "Build",
    name: "Storage footprint",
    metric: "total_size_bytes",
    color: "#cf5b56"
  }
];

const resourceAspects = [
  { key: "entity_count", label: "Entities", short: "Ent", color: "#4169a8" },
  { key: "association_count", label: "Associations", short: "Assoc", color: "#d49a2a" },
  { key: "interaction_count", label: "Interactions", short: "Inter", color: "#12a9b6" },
  { key: "identifier_count", label: "Identifiers", short: "IDs", color: "#7560a8" },
  { key: "ontology_term_count", label: "Ontology terms", short: "Onto", color: "#63a64e" }
];

function formatNumber(value) {
  return d3.format(",")(Number(value || 0));
}

function formatCompact(value) {
  return d3.format(".3~s")(Number(value || 0)).replace("G", "B");
}

function formatBytes(bytes) {
  const value = Number(bytes || 0);
  if (value >= 1024 ** 3) return `${d3.format(".2~f")(value / 1024 ** 3)} GB`;
  if (value >= 1024 ** 2) return `${d3.format(".2~f")(value / 1024 ** 2)} MB`;
  if (value >= 1024) return `${d3.format(".2~f")(value / 1024)} KB`;
  return `${formatNumber(value)} B`;
}

function formatMetric(value, metric) {
  return metric === "total_size_bytes" ? formatBytes(value) : formatNumber(value);
}

function termLabel(value) {
  if (value == null || value === "") return "Unspecified";
  const text = String(value);
  const match = text.match(/^[A-Z][A-Z0-9_]*:\d+:(.+)$/);
  if (match) return match[1];
  const omMatch = text.match(/^OM:\d+:(.+)$/);
  if (omMatch) return omMatch[1];
  return text;
}

function formatFacetLabel(value, facet) {
  if (facet === "taxonomy_id") {
    return taxonNames.has(String(value)) ? `${taxonNames.get(String(value))} (${value})` : `Taxon ${value}`;
  }
  if (facet === "source" || facet === "ontology_id") return String(value).replaceAll("_", " ");
  return termLabel(value);
}

function truncate(text, max = 32) {
  const string = String(text);
  return string.length > max ? `${string.slice(0, max - 1)}...` : string;
}

function clear(selector) {
  d3.select(selector).selectAll("*").remove();
}

function chartWidth(selector) {
  const node = d3.select(selector).node();
  return Math.max(320, Math.floor(node.getBoundingClientRect().width));
}

function tooltipHtml(title, rows) {
  const body = rows.map(([label, value]) => `<div>${label}: <strong>${value}</strong></div>`).join("");
  return `<strong>${title}</strong>${body}`;
}

const tooltip = d3.select("#tooltip");

function showTooltip(event, html) {
  tooltip
    .html(html)
    .style("left", `${event.clientX}px`)
    .style("top", `${event.clientY}px`)
    .style("opacity", 1);
}

function moveTooltip(event) {
  tooltip.style("left", `${event.clientX}px`).style("top", `${event.clientY}px`);
}

function hideTooltip() {
  tooltip.style("opacity", 0);
}

function renderKpis(data) {
  const items = [
    ["Entities", data.totals.entity_count, formatNumber],
    ["Relations", data.totals.relation_count, formatNumber],
    ["Interactions", data.totals.interaction_count, formatNumber],
    ["Associations", data.totals.association_count, formatNumber],
    ["Built resources", data.totals.built_resource_count, formatNumber],
    ["Storage", data.totals.database_size_bytes, formatBytes]
  ];

  d3.select("#generated-at").text(data.generated_at || "-");

  d3.select("#kpis")
    .selectAll(".kpi")
    .data(items)
    .join("div")
    .attr("class", "kpi")
    .html(([label, value, formatter]) => `<span>${label}</span><strong>${formatter(value)}</strong>`);
}

function resourceShortName(resource) {
  const name = resource.resource_name || resource.resource_id;
  return name
    .replace("Human Phenotype Ontology", "HPO")
    .replace("Mondo Disease Ontology", "Mondo")
    .replace("Guide to Pharmacology", "GtoP")
    .replace("ConnectomeDB2025", "Connectome")
    .replace("Phenol-Explorer", "PhenolExp");
}

function landscapeResourceColor(resource, row) {
  if (resource.resource_kind === "ontology") return palette.violet;
  if (row.metric === "interaction_count") return "#08a9b7";
  if (row.metric === "association_count") return "#d29226";
  if (row.metric === "ontology_term_count") return "#65a348";
  if (row.metric === "identifier_count") return "#7560a8";
  if (row.metric === "total_size_bytes") return "#c95d4f";
  return "#4169a8";
}

function proportionalSegments(items, width, minWidth = 4) {
  const total = d3.sum(items, (item) => item.sizeValue ?? item.value) || 1;
  const raw = items.map((item) => ((item.sizeValue ?? item.value) / total) * width);
  const adjusted = raw.map((itemWidth) => Math.max(minWidth, itemWidth));
  const scale = width / d3.sum(adjusted);
  let x = 0;

  return items.map((item, index) => {
    const segmentWidth = adjusted[index] * scale;
    const segment = { ...item, x, width: segmentWidth };
    x += segmentWidth;
    return segment;
  });
}

function drawLandscapeMosaic(svg, rows, resources, layout) {
  const {
    groupX,
    rowLabelX,
    plotX,
    top,
    plotWidth,
    rowHeight,
    rowGap
  } = layout;
  const groupRanges = [];

  rows.forEach((row, rowIndex) => {
    const y = top + rowIndex * (rowHeight + rowGap);
    const rowResources = resources
      .map((resource) => ({
        resource,
        value: Number(resource[row.metric] || 0)
      }))
      .filter((item) => item.value > 0)
      .map((item) => ({
        ...item,
        sizeValue: Math.log10(item.value + 1)
      }))
      .sort((a, b) => d3.descending(a.value, b.value) || d3.ascending(a.resource.resource_id, b.resource.resource_id));
    const segments = proportionalSegments(rowResources, plotWidth);

    svg.append("rect")
      .attr("x", plotX)
      .attr("y", y)
      .attr("width", plotWidth)
      .attr("height", rowHeight)
      .attr("fill", "#f1f4f2")
      .attr("stroke", "#dbe3df");

    svg.append("text")
      .attr("x", rowLabelX)
      .attr("y", y + rowHeight / 2)
      .attr("dy", "0.35em")
      .attr("text-anchor", "end")
      .attr("fill", "#17201c")
      .attr("font-size", 12)
      .attr("font-weight", 750)
      .text(row.name);

    const rowGroup = svg.append("g")
      .selectAll("g")
      .data(segments)
      .join("g")
      .attr("transform", (d) => `translate(${plotX + d.x},${y})`);

    rowGroup.append("rect")
      .attr("width", (d) => Math.max(0, d.width - 1))
      .attr("height", rowHeight)
      .attr("fill", (d) => landscapeResourceColor(d.resource, row))
      .attr("stroke", "#ffffff")
      .attr("stroke-width", 1)
      .on("mouseenter", (event, d) => showTooltip(event, tooltipHtml(resourceShortName(d.resource), [
        ["Layer", row.name],
        ["Value", formatMetric(d.value, row.metric)],
        ["Entities", formatNumber(d.resource.entity_count)],
        ["Relations", formatNumber(Number(d.resource.association_count || 0) + Number(d.resource.interaction_count || 0))]
      ])))
      .on("mousemove", moveTooltip)
      .on("mouseleave", hideTooltip);

    rowGroup.append("text")
      .attr("x", (d) => d.width / 2)
      .attr("y", rowHeight / 2)
      .attr("dy", "0.35em")
      .attr("text-anchor", "middle")
      .attr("fill", "#ffffff")
      .attr("font-size", (d) => d.width > 86 ? 10.5 : 9)
      .attr("font-weight", 750)
      .text((d) => {
        if (d.width < 34) return "";
        return truncate(resourceShortName(d.resource), Math.max(3, Math.floor(d.width / 7)));
      });

    const lastRange = groupRanges[groupRanges.length - 1];
    if (lastRange?.group === row.group) {
      lastRange.end = rowIndex;
    } else {
      groupRanges.push({ group: row.group, start: rowIndex, end: rowIndex });
    }
  });

  svg.selectAll(".landscape-group-label")
    .data(groupRanges)
    .join("text")
    .attr("class", "landscape-group-label")
    .attr("x", groupX)
    .attr("y", (d) => {
      const y0 = top + d.start * (rowHeight + rowGap);
      const y1 = top + d.end * (rowHeight + rowGap) + rowHeight;
      return (y0 + y1) / 2;
    })
    .attr("text-anchor", "middle")
    .attr("dominant-baseline", "middle")
    .attr("fill", "#66736e")
    .attr("font-size", 11)
    .attr("font-weight", 850)
    .attr("transform", (d) => {
      const y0 = top + d.start * (rowHeight + rowGap);
      const y1 = top + d.end * (rowHeight + rowGap) + rowHeight;
      const cy = (y0 + y1) / 2;
      return `rotate(-90 ${groupX} ${cy})`;
    })
    .text((d) => d.group.toUpperCase());
}

function miniPanel(svg, x, y, width, height, title) {
  const panel = svg.append("g").attr("transform", `translate(${x},${y})`);

  panel.append("text")
    .attr("x", 0)
    .attr("y", 0)
    .attr("fill", "#17201c")
    .attr("font-size", 12)
    .attr("font-weight", 850)
    .text(title);

  panel.append("line")
    .attr("x1", 0)
    .attr("x2", width)
    .attr("y1", 10)
    .attr("y2", 10)
    .attr("stroke", "#dce3df");

  return panel.append("g").attr("transform", "translate(0,20)");
}

function drawMiniBuildStatus(svg, resources, x, y, width, height) {
  const panel = miniPanel(svg, x, y, width, height, "A) Build status");
  const chartHeight = height - 38;
  const rows = [
    {
      label: "Data resources",
      success: resources.filter((r) => r.resource_kind !== "ontology" && r.build_status === "success").length,
      not_built: resources.filter((r) => r.resource_kind !== "ontology" && r.build_status !== "success").length
    },
    {
      label: "Ontologies",
      success: resources.filter((r) => r.resource_kind === "ontology" && r.build_status === "success").length,
      not_built: resources.filter((r) => r.resource_kind === "ontology" && r.build_status !== "success").length
    }
  ];
  const xScale = d3.scaleBand().domain(rows.map((d) => d.label)).range([36, width - 8]).padding(0.34);
  const yScale = d3.scaleLinear().domain([0, d3.max(rows, (d) => d.success + d.not_built) || 1]).nice().range([chartHeight, 8]);

  panel.append("g")
    .attr("transform", `translate(0,${chartHeight})`)
    .attr("class", "axis")
    .call(d3.axisBottom(xScale).tickSize(0))
    .call((axis) => axis.selectAll("text").attr("font-size", 10))
    .call((axis) => axis.select(".domain").remove());

  panel.append("g")
    .attr("class", "axis")
    .attr("transform", "translate(30,0)")
    .call(d3.axisLeft(yScale).ticks(3).tickSize(-width + 40))
    .call((axis) => axis.select(".domain").remove());

  const stacks = panel.selectAll(".build-stack")
    .data(rows)
    .join("g")
    .attr("transform", (d) => `translate(${xScale(d.label)},0)`);

  stacks.append("rect")
    .attr("y", (d) => yScale(d.success))
    .attr("width", xScale.bandwidth())
    .attr("height", (d) => chartHeight - yScale(d.success))
    .attr("fill", palette.green);

  stacks.append("rect")
    .attr("y", (d) => yScale(d.success + d.not_built))
    .attr("width", xScale.bandwidth())
    .attr("height", (d) => yScale(d.success) - yScale(d.success + d.not_built))
    .attr("fill", palette.coral);

  panel.append("text")
    .attr("x", width - 92)
    .attr("y", 8)
    .attr("font-size", 10)
    .attr("fill", palette.green)
    .attr("font-weight", 800)
    .text("built");

  panel.append("text")
    .attr("x", width - 52)
    .attr("y", 8)
    .attr("font-size", 10)
    .attr("fill", palette.coral)
    .attr("font-weight", 800)
    .text("pending");
}

function drawMiniLayerCoverage(svg, resources, x, y, width, height) {
  const panel = miniPanel(svg, x, y, width, height, "B) Resources by layer");
  const chartHeight = height - 38;
  const built = resources.filter((resource) => resource.build_status === "success");
  const layers = landscapeRows
    .filter((row) => row.metric !== "total_size_bytes")
    .map((row) => ({
      label: {
        association_count: "Assoc",
        interaction_count: "Inter",
        entity_count: "Entity",
        identifier_count: "IDs",
        ontology_term_count: "Onto"
      }[row.metric],
      value: built.filter((resource) => Number(resource[row.metric] || 0) > 0).length,
      color: row.color
    }));
  const xScale = d3.scaleBand().domain(layers.map((d) => d.label)).range([30, width - 4]).padding(0.26);
  const yScale = d3.scaleLinear().domain([0, d3.max(layers, (d) => d.value) || 1]).nice().range([chartHeight, 8]);

  panel.append("g")
    .attr("class", "axis")
    .attr("transform", `translate(0,${chartHeight})`)
    .call(d3.axisBottom(xScale).tickSize(0))
    .call((axis) => axis.selectAll("text").attr("font-size", 9))
    .call((axis) => axis.select(".domain").remove());

  panel.append("g")
    .attr("class", "axis")
    .attr("transform", "translate(26,0)")
    .call(d3.axisLeft(yScale).ticks(3).tickSize(-width + 34))
    .call((axis) => axis.select(".domain").remove());

  panel.selectAll("rect")
    .data(layers)
    .join("rect")
    .attr("x", (d) => xScale(d.label))
    .attr("y", (d) => yScale(d.value))
    .attr("width", xScale.bandwidth())
    .attr("height", (d) => chartHeight - yScale(d.value))
    .attr("rx", 3)
    .attr("fill", (d) => d.color);
}

function drawMiniRecordVolume(svg, data, x, y, width, height) {
  const panel = miniPanel(svg, x, y, width, height, "C) Records by layer");
  const chartHeight = height - 38;
  const totals = [
    ["Ent", data.totals.entity_count, "#4169a8"],
    ["Assoc", data.totals.association_count, "#d49a2a"],
    ["Inter", data.totals.interaction_count, "#12a9b6"],
    ["IDs", d3.sum(data.resources, (resource) => Number(resource.identifier_count || 0)), "#7560a8"],
    ["Onto", data.totals.ontology_term_count, "#63a64e"]
  ].map(([label, value, color]) => ({ label, value, color }));
  const xScale = d3.scaleBand().domain(totals.map((d) => d.label)).range([38, width - 4]).padding(0.28);
  const yScale = d3.scaleLog().domain([1, d3.max(totals, (d) => d.value) || 10]).range([chartHeight, 8]);

  panel.append("g")
    .attr("class", "axis")
    .attr("transform", `translate(0,${chartHeight})`)
    .call(d3.axisBottom(xScale).tickSize(0))
    .call((axis) => axis.selectAll("text").attr("font-size", 9))
    .call((axis) => axis.select(".domain").remove());

  panel.append("g")
    .attr("class", "axis")
    .attr("transform", "translate(32,0)")
    .call(d3.axisLeft(yScale).ticks(4).tickFormat((d) => formatCompact(d)).tickSize(-width + 40))
    .call((axis) => axis.select(".domain").remove());

  panel.selectAll("rect")
    .data(totals)
    .join("rect")
    .attr("x", (d) => xScale(d.label))
    .attr("y", (d) => yScale(d.value))
    .attr("width", xScale.bandwidth())
    .attr("height", (d) => chartHeight - yScale(d.value))
    .attr("rx", 3)
    .attr("fill", (d) => d.color)
    .on("mouseenter", (event, d) => showTooltip(event, tooltipHtml(d.label, [["Records", formatNumber(d.value)]])))
    .on("mousemove", moveTooltip)
    .on("mouseleave", hideTooltip);
}

function drawMiniLayerOverlap(svg, resources, x, y, width, height) {
  const panel = miniPanel(svg, x, y, width, height, "D) Resource overlap");
  const chartHeight = height - 38;
  const contentMetrics = landscapeRows
    .filter((row) => row.metric !== "total_size_bytes")
    .map((row) => row.metric);
  const overlap = d3.rollups(
    resources.filter((resource) => resource.build_status === "success"),
    (items) => items.length,
    (resource) => contentMetrics.filter((metric) => Number(resource[metric] || 0) > 0).length
  )
    .map(([layers, value]) => ({ layers, value }))
    .sort((a, b) => d3.ascending(a.layers, b.layers));
  const xScale = d3.scaleBand().domain(overlap.map((d) => d.layers)).range([32, width - 6]).padding(0.28);
  const yScale = d3.scaleLinear().domain([0, d3.max(overlap, (d) => d.value) || 1]).nice().range([chartHeight, 8]);
  const color = d3.scaleOrdinal().domain(overlap.map((d) => d.layers)).range(["#4169a8", "#12a9b6", "#63a64e", "#d49a2a", "#c95d4f"]);

  panel.append("g")
    .attr("class", "axis")
    .attr("transform", `translate(0,${chartHeight})`)
    .call(d3.axisBottom(xScale).tickSize(0).tickFormat((d) => `${d}L`))
    .call((axis) => axis.selectAll("text").attr("font-size", 9))
    .call((axis) => axis.select(".domain").remove());

  panel.append("g")
    .attr("class", "axis")
    .attr("transform", "translate(26,0)")
    .call(d3.axisLeft(yScale).ticks(3).tickSize(-width + 34))
    .call((axis) => axis.select(".domain").remove());

  panel.selectAll("rect")
    .data(overlap)
    .join("rect")
    .attr("x", (d) => xScale(d.layers))
    .attr("y", (d) => yScale(d.value))
    .attr("width", xScale.bandwidth())
    .attr("height", (d) => chartHeight - yScale(d.value))
    .attr("rx", 3)
    .attr("fill", (d) => color(d.layers));
}

function aspectMaxima(resources) {
  return Object.fromEntries(
    resourceAspects.map((aspect) => [
      aspect.key,
      d3.max(resources, (resource) => Number(resource[aspect.key] || 0)) || 1
    ])
  );
}

function normalizedAspectValue(resource, aspect, maxima) {
  const value = Number(resource[aspect.key] || 0);
  const max = maxima[aspect.key] || 1;
  return Math.log10(value + 1) / Math.log10(max + 1);
}

function resourceFingerprintScore(resource, maxima) {
  return d3.sum(resourceAspects, (aspect) => normalizedAspectValue(resource, aspect, maxima));
}

function dominantResourceAspect(resource, maxima) {
  return resourceAspects
    .filter((aspect) => aspect.key !== "total_size_bytes")
    .map((aspect) => ({
      ...aspect,
      normalized: normalizedAspectValue(resource, aspect, maxima)
    }))
    .sort((a, b) => d3.descending(a.normalized, b.normalized))[0];
}

function totalContentRecords(resource) {
  return d3.sum(resourceAspects.filter((aspect) => aspect.key !== "total_size_bytes"), (aspect) => Number(resource[aspect.key] || 0));
}

function aggregateResourceGroup(resources, label) {
  const aggregate = {
    resource_id: "other_sources",
    resource_name: label,
    resource_kind: "aggregate",
    build_status: "success",
    is_aggregate: true,
    source_count: resources.length
  };

  for (const aspect of resourceAspects) {
    aggregate[aspect.key] = d3.sum(resources, (resource) => Number(resource[aspect.key] || 0));
  }

  aggregate.association_count = d3.sum(resources, (resource) => Number(resource.association_count || 0));
  aggregate.interaction_count = d3.sum(resources, (resource) => Number(resource.interaction_count || 0));
  aggregate.entity_count = d3.sum(resources, (resource) => Number(resource.entity_count || 0));
  aggregate.identifier_count = d3.sum(resources, (resource) => Number(resource.identifier_count || 0));
  aggregate.ontology_term_count = d3.sum(resources, (resource) => Number(resource.ontology_term_count || 0));
  aggregate.total_size_bytes = d3.sum(resources, (resource) => Number(resource.total_size_bytes || 0));

  return aggregate;
}

function makeResourceFingerprintData(resources, topLimit) {
  const maxima = aspectMaxima(resources);
  const scored = resources
    .map((resource) => {
      const dominant = dominantResourceAspect(resource, maxima);
      return {
        ...resource,
        score: resourceFingerprintScore(resource, maxima),
        dominant_aspect: dominant?.label || "Resources",
        dominant_color: dominant?.color || palette.gray
      };
    })
    .sort((a, b) => d3.descending(a.score, b.score) || d3.ascending(a.resource_id, b.resource_id));
  const top = scored.slice(0, topLimit);
  const rest = scored.slice(topLimit);
  const tiles = rest.length ? [...top, aggregateResourceGroup(rest, `Other built sources (${rest.length})`)] : top;
  const tileMaxima = aspectMaxima(tiles);

  return {
    maxima: tileMaxima,
    tiles: tiles.map((resource) => {
      const dominant = resource.is_aggregate
        ? { label: "Other sources", color: "#8a9690" }
        : dominantResourceAspect(resource, tileMaxima);
      return {
        ...resource,
        score: resourceFingerprintScore(resource, tileMaxima),
        dominant_aspect: dominant?.label || "Resources",
        dominant_color: dominant?.color || palette.gray
      };
    })
  };
}

function drawResourceFingerprintTile(svg, leaf, maxima) {
  const resource = leaf.data;
  const x = leaf.x0;
  const y = leaf.y0;
  const width = leaf.x1 - leaf.x0;
  const height = leaf.y1 - leaf.y0;
  const usable = width > 42 && height > 28;
  const clipId = `clip-${resource.resource_id.replace(/[^a-zA-Z0-9_-]/g, "-")}`;

  svg.append("clipPath")
    .attr("id", clipId)
    .append("rect")
    .attr("x", x + 5)
    .attr("y", y + 5)
    .attr("width", Math.max(0, width - 10))
    .attr("height", Math.max(0, height - 10))
    .attr("rx", 4);

  const tile = svg.append("g")
    .attr("clip-path", `url(#${clipId})`);

  tile.append("rect")
    .attr("x", x)
    .attr("y", y)
    .attr("width", width)
    .attr("height", height)
    .attr("rx", 5)
    .attr("fill", resource.dominant_color)
    .attr("opacity", resource.is_aggregate ? 0.82 : 0.94)
    .attr("stroke", "#ffffff")
    .attr("stroke-width", 1.4)
    .on("mouseenter", (event) => showTooltip(event, tooltipHtml(resource.resource_name, [
      ["Dominant aspect", resource.dominant_aspect],
      ["Entities", formatNumber(resource.entity_count)],
      ["Associations", formatNumber(resource.association_count)],
      ["Interactions", formatNumber(resource.interaction_count)],
      ["Identifiers", formatNumber(resource.identifier_count)],
      ["Ontology terms", formatNumber(resource.ontology_term_count)],
      ["Storage", formatBytes(resource.total_size_bytes)]
    ])))
    .on("mousemove", moveTooltip)
    .on("mouseleave", hideTooltip);

  if (!usable) return;

  const name = resourceShortName(resource);
  const contentTotal = totalContentRecords(resource);
  const fontSize = Math.max(10, Math.min(15, width / 12, height / 5));
  const stripHeight = Math.min(30, Math.max(18, height * 0.28));
  const labelY = y + Math.min(24, height * 0.28);

  tile.append("text")
    .attr("x", x + 8)
    .attr("y", labelY)
    .attr("fill", "#ffffff")
    .attr("font-size", fontSize)
    .attr("font-weight", 900)
    .text(truncate(name, Math.max(4, Math.floor((width - 14) / (fontSize * 0.58)))));

  if (height > 56 && width > 95) {
    tile.append("text")
      .attr("x", x + 8)
      .attr("y", labelY + fontSize + 7)
      .attr("fill", "rgba(255,255,255,0.86)")
      .attr("font-size", Math.max(9, fontSize - 3))
      .attr("font-weight", 750)
      .text(resource.is_aggregate ? `${resource.source_count} resources` : `${formatCompact(contentTotal)} records`);
  }

  if (height <= 48 || width <= 88) return;

  const stripY = y + height - stripHeight - 7;
  const stripX = x + 8;
  const stripWidth = width - 16;
  const slotWidth = stripWidth / resourceAspects.length;
  const maxBarHeight = stripHeight - 8;

  const strip = tile.append("g");
  resourceAspects.forEach((aspect, index) => {
    const normalized = normalizedAspectValue(resource, aspect, maxima);
    const barHeight = Math.max(2, normalized * maxBarHeight);
    const barWidth = Math.max(4, Math.min(14, slotWidth * 0.58));
    const barX = stripX + index * slotWidth + (slotWidth - barWidth) / 2;
    const barY = stripY + maxBarHeight - barHeight;

    strip.append("rect")
      .attr("x", barX)
      .attr("y", stripY)
      .attr("width", barWidth)
      .attr("height", maxBarHeight)
      .attr("rx", 2)
      .attr("fill", "rgba(255,255,255,0.22)");

    strip.append("rect")
      .attr("x", barX)
      .attr("y", barY)
      .attr("width", barWidth)
      .attr("height", barHeight)
      .attr("rx", 2)
      .attr("fill", aspect.color);
  });
}

function drawResourceFingerprintLegend(svg, x, y, width) {
  svg.append("text")
    .attr("x", x)
    .attr("y", y)
    .attr("fill", "#17201c")
    .attr("font-size", 12)
    .attr("font-weight", 850)
    .text("Tile area combines all aspects; bars inside each tile show log-normalized aspect strength");

  const columns = width < 760 ? 3 : 6;
  const itemWidth = width / columns;
  const legend = svg.append("g").attr("transform", `translate(${x},${y + 18})`);

  resourceAspects.forEach((aspect, index) => {
    const item = legend.append("g")
      .attr("transform", `translate(${(index % columns) * itemWidth},${Math.floor(index / columns) * 20})`);

    item.append("rect")
      .attr("x", 0)
      .attr("y", -9)
      .attr("width", 11)
      .attr("height", 11)
      .attr("rx", 2)
      .attr("fill", aspect.color);

    item.append("text")
      .attr("x", 17)
      .attr("y", 0)
      .attr("fill", "#66736e")
      .attr("font-size", 11)
      .attr("font-weight", 750)
      .text(aspect.label);
  });
}

function stableResourceColorScale(resources) {
  const domain = resources
    .map((resource) => resource.resource_id)
    .sort(d3.ascending);
  const colors = domain.map((_, index) => d3.hsl((index * 137.508) % 360, 0.48, 0.48).formatHex());

  return d3.scaleOrdinal(domain, colors);
}

function drawAspectResourceTreemap(panel, aspect, resources, colorScale, width, height) {
  const headerHeight = 28;
  const aspectResources = resources
    .map((resource) => ({
      ...resource,
      aspect_value: Number(resource[aspect.key] || 0),
      log_value: Math.log10(Number(resource[aspect.key] || 0) + 1)
    }))
    .filter((resource) => resource.aspect_value > 0)
    .sort((a, b) => d3.descending(a.aspect_value, b.aspect_value) || d3.ascending(a.resource_id, b.resource_id));

  panel.append("rect")
    .attr("x", 0)
    .attr("y", 0)
    .attr("width", width)
    .attr("height", height)
    .attr("rx", 7)
    .attr("fill", "#f2f5f3")
    .attr("stroke", "#dce3df");

  panel.append("text")
    .attr("x", 10)
    .attr("y", 18)
    .attr("fill", aspect.color)
    .attr("font-size", 12)
    .attr("font-weight", 950)
    .text(aspect.label.toUpperCase());

  panel.append("text")
    .attr("x", width - 10)
    .attr("y", 18)
    .attr("text-anchor", "end")
    .attr("fill", "#66736e")
    .attr("font-size", 10.5)
    .attr("font-weight", 750)
    .text(`${aspectResources.length} resources`);

  if (!aspectResources.length) return;

  const root = d3.hierarchy({ children: aspectResources })
    .sum((resource) => resource.log_value || 0)
    .sort((a, b) => d3.descending(a.value, b.value));

  d3.treemap()
    .size([width - 12, height - headerHeight - 10])
    .paddingInner(2)
    .round(true)(root);

  const plot = panel.append("g").attr("transform", `translate(6,${headerHeight})`);

  const tiles = plot.selectAll("g")
    .data(root.leaves())
    .join("g")
    .attr("transform", (leaf) => `translate(${leaf.x0},${leaf.y0})`);

  tiles.append("rect")
    .attr("width", (leaf) => Math.max(0, leaf.x1 - leaf.x0))
    .attr("height", (leaf) => Math.max(0, leaf.y1 - leaf.y0))
    .attr("rx", 4)
    .attr("fill", (leaf) => colorScale(leaf.data.resource_id))
    .attr("stroke", "#ffffff")
    .attr("stroke-width", 1)
    .on("mouseenter", (event, leaf) => showTooltip(event, tooltipHtml(resourceShortName(leaf.data), [
      [aspect.label, formatMetric(leaf.data.aspect_value, aspect.key)],
      ["Entities", formatNumber(leaf.data.entity_count)],
      ["Associations", formatNumber(leaf.data.association_count)],
      ["Interactions", formatNumber(leaf.data.interaction_count)],
      ["Identifiers", formatNumber(leaf.data.identifier_count)],
      ["Ontology terms", formatNumber(leaf.data.ontology_term_count)],
      ["Storage", formatBytes(leaf.data.total_size_bytes)]
    ])))
    .on("mousemove", moveTooltip)
    .on("mouseleave", hideTooltip);

  tiles.each(function (leaf) {
    const tile = d3.select(this);
    const tileWidth = leaf.x1 - leaf.x0;
    const tileHeight = leaf.y1 - leaf.y0;
    if (tileWidth < 28 || tileHeight < 18) return;

    const name = resourceShortName(leaf.data);
    const fontSize = Math.max(8, Math.min(12, tileWidth / 8, tileHeight / 3.2));
    const maxChars = Math.max(2, Math.floor((tileWidth - 8) / (fontSize * 0.58)));

    tile.append("text")
      .attr("x", 5)
      .attr("y", Math.min(15, tileHeight - 4))
      .attr("fill", "#ffffff")
      .attr("font-size", fontSize)
      .attr("font-weight", 900)
      .text(truncate(name, maxChars));

    if (tileWidth >= 74 && tileHeight >= 42) {
      tile.append("text")
        .attr("x", 5)
        .attr("y", Math.min(30, tileHeight - 7))
        .attr("fill", "rgba(255,255,255,0.86)")
        .attr("font-size", Math.max(8, fontSize - 2))
        .attr("font-weight", 750)
        .text(formatMetric(leaf.data.aspect_value, aspect.key));
    }
  });
}

function drawResourceColorLegend(svg, resources, colorScale, x, y, width) {
  const visible = resources
    .slice()
    .sort((a, b) => d3.descending(totalContentRecords(a), totalContentRecords(b)))
    .slice(0, width < 760 ? 12 : 18);
  const columns = width < 760 ? 3 : 6;
  const itemWidth = width / columns;

  svg.append("text")
    .attr("x", x)
    .attr("y", y)
    .attr("fill", "#17201c")
    .attr("font-size", 12)
    .attr("font-weight", 850)
    .text("Same resource, same color across every aspect");

  const legend = svg.append("g").attr("transform", `translate(${x},${y + 20})`);

  visible.forEach((resource, index) => {
    const item = legend.append("g")
      .attr("transform", `translate(${(index % columns) * itemWidth},${Math.floor(index / columns) * 18})`);

    item.append("rect")
      .attr("x", 0)
      .attr("y", -9)
      .attr("width", 10)
      .attr("height", 10)
      .attr("rx", 2)
      .attr("fill", colorScale(resource.resource_id));

    item.append("text")
      .attr("x", 16)
      .attr("y", 0)
      .attr("fill", "#66736e")
      .attr("font-size", 10.5)
      .attr("font-weight", 750)
      .text(truncate(resourceShortName(resource), Math.max(7, Math.floor(itemWidth / 8))));
  });
}

function renderResourceLandscape(data) {
  clear("#resource-landscape");

  const resources = data.resources.filter((resource) => resource.build_status === "success");
  if (!resources.length) {
    d3.select("#resource-landscape").append("div").attr("class", "empty-state").text("No built resources");
    return;
  }

  const containerWidth = chartWidth("#resource-landscape");
  const compactFigure = containerWidth < 860;
  const width = Math.max(680, containerWidth);
  const margin = { top: 68, right: 16, bottom: 18, left: 16 };
  const plotWidth = width - margin.left - margin.right;
  const panelGap = 8;
  const lowerColumns = compactFigure ? 2 : 4;
  const lowerAspects = resourceAspects.filter((aspect) => aspect.key !== "entity_count");
  const entityHeight = compactFigure ? 260 : 250;
  const lowerPanelWidth = (plotWidth - panelGap * (lowerColumns - 1)) / lowerColumns;
  const lowerPanelHeight = compactFigure ? 220 : 230;
  const lowerRows = Math.ceil(lowerAspects.length / lowerColumns);
  const plotHeight = entityHeight + panelGap + lowerRows * lowerPanelHeight + (lowerRows - 1) * panelGap;
  const legendHeight = compactFigure ? 92 : 72;
  const height = margin.top + plotHeight + legendHeight + margin.bottom;
  const colorScale = stableResourceColorScale(resources);

  const svg = d3.select("#resource-landscape")
    .append("svg")
    .attr("viewBox", `0 0 ${width} ${height}`)
    .attr("role", "img")
    .attr("aria-label", "OmniPath database resource fingerprints");

  svg.append("text")
    .attr("x", width / 2)
    .attr("y", 22)
    .attr("text-anchor", "middle")
    .attr("fill", "#17201c")
    .attr("font-size", 18)
    .attr("font-weight", 900)
    .text("OmniPath Database Resources by Aspect");

  svg.append("text")
    .attr("x", width / 2)
    .attr("y", 42)
    .attr("text-anchor", "middle")
    .attr("fill", "#66736e")
    .attr("font-size", 11.5)
    .attr("font-weight", 650)
    .text("Resource tiles are sized by log10(value + 1) and keep a stable color across aspects");

  const plot = svg.append("g").attr("transform", `translate(${margin.left},${margin.top})`);
  const entityAspect = resourceAspects.find((aspect) => aspect.key === "entity_count");

  drawAspectResourceTreemap(plot.append("g"), entityAspect, resources, colorScale, plotWidth, entityHeight);

  lowerAspects.forEach((aspect, index) => {
    const column = index % lowerColumns;
    const row = Math.floor(index / lowerColumns);
    const panel = plot.append("g")
      .attr("transform", `translate(${column * (lowerPanelWidth + panelGap)},${entityHeight + panelGap + row * (lowerPanelHeight + panelGap)})`);

    drawAspectResourceTreemap(panel, aspect, resources, colorScale, lowerPanelWidth, lowerPanelHeight);
  });

  drawResourceColorLegend(svg, resources, colorScale, margin.left, margin.top + plotHeight + 28, plotWidth);
}

function drawHorizontalBars(selector, rows, options = {}) {
  clear(selector);

  const data = rows.filter((row) => Number(row.value) > 0);
  if (!data.length) {
    d3.select(selector).append("div").attr("class", "empty-state").text("No data");
    return;
  }

  const width = chartWidth(selector);
  const compact = width < 620;
  const margin = {
    top: 12,
    right: compact ? 70 : 92,
    bottom: 30,
    left: compact ? 128 : 210
  };
  const rowHeight = compact ? 24 : 28;
  const height = margin.top + margin.bottom + data.length * rowHeight;

  const svg = d3.select(selector).append("svg").attr("viewBox", `0 0 ${width} ${height}`);
  const innerWidth = Math.max(60, width - margin.left - margin.right);
  const innerHeight = height - margin.top - margin.bottom;
  const g = svg.append("g").attr("transform", `translate(${margin.left},${margin.top})`);

  const x = d3.scaleLinear().domain([0, d3.max(data, (d) => d.value) || 1]).nice().range([0, innerWidth]);
  const y = d3.scaleBand().domain(data.map((d) => d.key)).range([0, innerHeight]).padding(0.24);

  g.append("g")
    .attr("class", "grid")
    .attr("transform", `translate(0,${innerHeight})`)
    .call(d3.axisBottom(x).ticks(compact ? 3 : 5).tickSize(-innerHeight).tickFormat(""))
    .call((axis) => axis.select(".domain").remove());

  g.selectAll("rect")
    .data(data)
    .join("rect")
    .attr("x", 0)
    .attr("y", (d) => y(d.key))
    .attr("width", (d) => x(d.value))
    .attr("height", y.bandwidth())
    .attr("rx", 3)
    .attr("fill", (d) => d.color || palette.teal)
    .on("mouseenter", (event, d) => showTooltip(event, d.tooltip))
    .on("mousemove", moveTooltip)
    .on("mouseleave", hideTooltip);

  g.append("g")
    .attr("class", "axis")
    .attr("transform", `translate(0,${innerHeight})`)
    .call(d3.axisBottom(x).ticks(compact ? 3 : 5).tickFormat((d) => formatCompact(d)))
    .call((axis) => axis.select(".domain").remove());

  g.selectAll(".bar-label")
    .data(data)
    .join("text")
    .attr("class", "bar-label")
    .attr("x", -10)
    .attr("y", (d) => y(d.key) + y.bandwidth() / 2)
    .attr("dy", "0.35em")
    .attr("text-anchor", "end")
    .text((d) => truncate(d.label, compact ? 17 : 30))
    .append("title")
    .text((d) => d.label);

  g.selectAll(".value-label")
    .data(data)
    .join("text")
    .attr("class", "value-label")
    .attr("x", (d) => Math.min(innerWidth + 8, x(d.value) + 7))
    .attr("y", (d) => y(d.key) + y.bandwidth() / 2)
    .attr("dy", "0.35em")
    .text((d) => options.valueFormatter ? options.valueFormatter(d.value, d) : formatCompact(d.value));
}

function renderResourceBars(data) {
  const metric = state.resourceMetric;
  const rows = data.resources
    .filter((resource) => resource.build_status === "success")
    .map((resource) => ({
      key: resource.resource_id,
      label: resource.resource_name,
      value: Number(resource[metric] || 0),
      color: resource.resource_kind === "ontology" ? palette.violet : palette.teal,
      tooltip: tooltipHtml(resource.resource_name, [
        ["Metric", metricLabels[metric]],
        ["Value", formatMetric(resource[metric], metric)],
        ["Kind", resource.resource_kind || "unknown"],
        ["Built", resource.last_built_at ? new Date(resource.last_built_at).toLocaleString() : "-"]
      ])
    }))
    .sort((a, b) => d3.descending(a.value, b.value))
    .slice(0, 20);

  drawHorizontalBars("#resource-bars", rows, {
    valueFormatter: (value) => formatMetric(value, metric)
  });
}

function sankeyNodeLabel(id) {
  return id.replace(/^(subject|predicate|object):/, "");
}

function sankeyNodeLayer(id) {
  if (id.startsWith("subject:")) return "Subject type";
  if (id.startsWith("predicate:")) return "Predicate";
  return "Object type";
}

function sankeyNodeColor(node) {
  if (node.id.startsWith("subject:")) return palette.teal;
  if (node.id.startsWith("object:")) return palette.violet;
  return node.category === "interaction" ? palette.blue : node.category === "association" ? palette.amber : palette.gray;
}

function buildRelationFlowGraph(rows) {
  const nodeMap = new Map();
  const linkMap = new Map();

  function ensureNode(id, category = null) {
    if (!nodeMap.has(id)) {
      nodeMap.set(id, { id, name: sankeyNodeLabel(id), category, relationCount: 0 });
    } else if (category && !nodeMap.get(id).category) {
      nodeMap.get(id).category = category;
    }
  }

  function addLink(source, target, value, weight, category) {
    const key = `${source}->${target}`;
    const existing = linkMap.get(key);
    if (existing) {
      existing.value += weight;
      existing.relationCount += value;
      existing.category = existing.category || category;
      return;
    }
    linkMap.set(key, { source, target, value: weight, relationCount: value, category });
  }

  for (const row of rows) {
    const value = Number(row.relation_count || 0);
    if (value <= 0) continue;
    const weight = Math.log10(value + 1);

    const subject = `subject:${formatFacetLabel(row.subject_type, "entity_type")}`;
    const predicate = `predicate:${formatFacetLabel(row.predicate, "predicate")}`;
    const object = `object:${formatFacetLabel(row.object_type, "entity_type")}`;

    ensureNode(subject);
    ensureNode(predicate, row.category);
    ensureNode(object);
    nodeMap.get(subject).relationCount += value;
    nodeMap.get(predicate).relationCount += value;
    nodeMap.get(object).relationCount += value;
    addLink(subject, predicate, value, weight, row.category);
    addLink(predicate, object, value, weight, row.category);
  }

  const nodes = Array.from(nodeMap.values());
  const nodeIds = new Set(nodes.map((node) => node.id));
  const links = Array.from(linkMap.values()).filter((link) => nodeIds.has(link.source) && nodeIds.has(link.target));

  return { nodes, links };
}

function renderRelationFlow(data) {
  clear("#relation-flow");

  const rows = (data.subject_predicate_object_types || [])
    .filter((row) => Number(row.relation_count || 0) > 0)
    .slice(0, 70);

  if (!rows.length || !d3.sankey) {
    d3.select("#relation-flow").append("div").attr("class", "empty-state").text("No relation flow data");
    return;
  }

  const graph = buildRelationFlowGraph(rows);
  const width = chartWidth("#relation-flow");
  const compact = width < 560;
  const nodeWidth = compact ? 14 : 18;
  const margin = { top: 36, right: compact ? 78 : 128, bottom: 20, left: compact ? 78 : 128 };
  const height = compact ? 720 : 820;
  const innerWidth = width - margin.left - margin.right;
  const innerHeight = height - margin.top - margin.bottom;
  const svg = d3.select("#relation-flow").append("svg").attr("viewBox", `0 0 ${width} ${height}`);

  svg.append("text")
    .attr("x", margin.left)
    .attr("y", 18)
    .attr("fill", "#17201c")
    .attr("font-size", 12)
    .attr("font-weight", 850)
    .text(`${formatNumber(d3.sum(rows, (row) => Number(row.relation_count || 0)))} relations across ${rows.length} top paths, widths scaled as log10(count + 1)`);

  const sankey = d3.sankey()
    .nodeId((node) => node.id)
    .nodeWidth(nodeWidth)
    .nodePadding(compact ? 7 : 10)
    .nodeAlign(d3.sankeyJustify)
    .extent([[margin.left, margin.top], [margin.left + innerWidth, margin.top + innerHeight]]);

  const layout = sankey({
    nodes: graph.nodes.map((node) => ({ ...node })),
    links: graph.links.map((link) => ({ ...link }))
  });

  const linkColor = (link) => link.category === "interaction" ? palette.blue : link.category === "association" ? palette.amber : palette.gray;

  svg.append("g")
    .attr("fill", "none")
    .selectAll("path")
    .data(layout.links)
    .join("path")
    .attr("d", d3.sankeyLinkHorizontal())
    .attr("stroke", linkColor)
    .attr("stroke-opacity", 0.34)
    .attr("stroke-width", (link) => Math.max(1, link.width))
    .on("mouseenter", (event, link) => {
      d3.select(event.currentTarget).attr("stroke-opacity", 0.68);
      showTooltip(event, tooltipHtml(`${link.source.name} -> ${link.target.name}`, [
        ["Relations", formatNumber(link.relationCount)],
        ["Scaled width", d3.format(".2f")(link.value)],
        ["Category", link.category || "-"]
      ]));
    })
    .on("mousemove", moveTooltip)
    .on("mouseleave", (event) => {
      d3.select(event.currentTarget).attr("stroke-opacity", 0.34);
      hideTooltip();
    });

  const node = svg.append("g")
    .selectAll("g")
    .data(layout.nodes)
    .join("g");

  node.append("rect")
    .attr("x", (d) => d.x0)
    .attr("y", (d) => d.y0)
    .attr("height", (d) => Math.max(1, d.y1 - d.y0))
    .attr("width", (d) => d.x1 - d.x0)
    .attr("rx", 2)
    .attr("fill", sankeyNodeColor)
    .attr("stroke", "#ffffff")
    .attr("stroke-width", 1)
    .on("mouseenter", (event, d) => showTooltip(event, tooltipHtml(d.name, [
      ["Layer", sankeyNodeLayer(d.id)],
      ["Relations", formatNumber(d.relationCount)],
      ["Scaled width", d3.format(".2f")(d.value)],
      ["Category", d.category || "-"]
    ])))
    .on("mousemove", moveTooltip)
    .on("mouseleave", hideTooltip);

  node.append("text")
    .attr("x", (d) => {
      if (d.id.startsWith("subject:")) return d.x0 - 8;
      if (d.id.startsWith("object:")) return d.x1 + 8;
      return (d.x0 + d.x1) / 2;
    })
    .attr("y", (d) => (d.y0 + d.y1) / 2)
    .attr("dy", "0.35em")
    .attr("text-anchor", (d) => {
      if (d.id.startsWith("subject:")) return "end";
      if (d.id.startsWith("object:")) return "start";
      return "middle";
    })
    .attr("class", "bar-label")
    .attr("font-size", compact ? 10 : 12)
    .style("fill", "#17201c")
    .style("paint-order", "stroke")
    .style("stroke", "rgba(255,255,255,0.88)")
    .style("stroke-width", 4)
    .style("stroke-linejoin", "round")
    .text((d) => {
      const nodeHeight = d.y1 - d.y0;
      if (nodeHeight < (compact ? 18 : 14)) return "";
      return truncate(d.name, d.id.startsWith("predicate:") ? 20 : compact ? 12 : 18);
    })
    .append("title")
    .text((d) => d.name);

  const headers = [
    ["Subject type", margin.left],
    ["Predicate", width / 2],
    ["Object type", width - margin.right]
  ];

  svg.append("g")
    .selectAll("text")
    .data(headers)
    .join("text")
    .attr("x", ([, x]) => x)
    .attr("y", 34)
    .attr("text-anchor", (_, index) => index === 0 ? "start" : index === 1 ? "middle" : "end")
    .attr("fill", "#66736e")
    .attr("font-size", 11)
    .attr("font-weight", 850)
    .text(([label]) => label.toUpperCase());
}

function renderEntityFacet(data) {
  const facet = state.entityFacet;
  const rows = (data.entity_facets[facet] || [])
    .slice(0, facet === "taxonomy_id" ? 24 : 22)
    .map((item) => ({
      key: item.label,
      label: formatFacetLabel(item.label, facet),
      value: Number(item.count || 0),
      color: categorical(facet),
      tooltip: tooltipHtml(formatFacetLabel(item.label, facet), [
        ["Facet", facet.replaceAll("_", " ")],
        ["Entities", formatNumber(item.count)]
      ])
    }));

  drawHorizontalBars("#entity-facets", rows);
}

function renderRelationFacet(data) {
  const facet = state.relationFacet;
  const rows = (data.relation_facets[facet] || [])
    .slice(0, 22)
    .map((item) => ({
      key: item.label,
      label: formatFacetLabel(item.label, facet),
      value: Number(item.count || 0),
      color:
        item.category === "interaction"
          ? palette.blue
          : item.category === "association"
            ? palette.amber
            : categorical(item.label),
      tooltip: tooltipHtml(formatFacetLabel(item.label, facet), [
        ["Facet", facet.replaceAll("_", " ")],
        ["Relations", formatNumber(item.count)],
        ["Category", item.category || "-"]
      ])
    }));

  drawHorizontalBars("#relation-facets", rows);
}

function renderResourceMatrix(data) {
  clear("#resource-matrix");

  const columns = [
    ["entity_count", "Entities"],
    ["association_count", "Associations"],
    ["interaction_count", "Interactions"],
    ["identifier_count", "Identifiers"],
    ["ontology_term_count", "Ontology terms"]
  ];

  const rows = data.resources
    .filter((resource) => resource.build_status === "success")
    .map((resource) => ({
      ...resource,
      profileTotal: columns.reduce((sum, [key]) => sum + Number(resource[key] || 0), 0)
    }))
    .sort((a, b) => d3.descending(a.profileTotal, b.profileTotal))
    .slice(0, 24);

  const width = chartWidth("#resource-matrix");
  const compact = width < 720;
  const margin = {
    top: compact ? 66 : 52,
    right: 18,
    bottom: 22,
    left: compact ? 104 : 174
  };
  const cellHeight = compact ? 23 : 26;
  const height = margin.top + margin.bottom + rows.length * cellHeight;
  const innerWidth = width - margin.left - margin.right;
  const cellWidth = innerWidth / columns.length;
  const values = rows.flatMap((row) => columns.map(([key]) => Number(row[key] || 0))).filter((value) => value > 0);
  const color = d3.scaleSequentialLog(d3.interpolateYlGnBu).domain([1, d3.max(values) || 10]);

  const svg = d3.select("#resource-matrix").append("svg").attr("viewBox", `0 0 ${width} ${height}`);
  const g = svg.append("g").attr("transform", `translate(${margin.left},${margin.top})`);

  svg.append("g")
    .attr("transform", `translate(${margin.left},${compact ? 12 : 18})`)
    .selectAll("text")
    .data(columns)
    .join("text")
    .attr("x", (_, index) => index * cellWidth + cellWidth / 2)
    .attr("y", 0)
    .attr("text-anchor", "middle")
    .attr("class", "bar-label")
    .text(([, label]) => compact ? truncate(label, 10) : label);

  svg.append("g")
    .attr("transform", `translate(${margin.left - 10},${margin.top})`)
    .selectAll("text")
    .data(rows)
    .join("text")
    .attr("x", 0)
    .attr("y", (_, index) => index * cellHeight + cellHeight / 2)
    .attr("dy", "0.35em")
    .attr("text-anchor", "end")
    .attr("class", "bar-label")
    .text((row) => truncate(row.resource_name, compact ? 13 : 24))
    .append("title")
    .text((row) => row.resource_name);

  const cell = g.selectAll("g.cell")
    .data(rows.flatMap((row) => columns.map(([key, label], columnIndex) => ({ row, key, label, columnIndex }))))
    .join("g")
    .attr("class", "cell")
    .attr("transform", (d, index) => {
      const rowIndex = Math.floor(index / columns.length);
      return `translate(${d.columnIndex * cellWidth},${rowIndex * cellHeight})`;
    });

  cell.append("rect")
    .attr("x", 1)
    .attr("y", 1)
    .attr("width", Math.max(1, cellWidth - 2))
    .attr("height", cellHeight - 2)
    .attr("rx", 3)
    .attr("fill", (d) => Number(d.row[d.key] || 0) > 0 ? color(Number(d.row[d.key])) : "#eef2f0")
    .on("mouseenter", (event, d) => showTooltip(event, tooltipHtml(d.row.resource_name, [
      [d.label, formatNumber(d.row[d.key])],
      ["Storage", formatBytes(d.row.total_size_bytes)]
    ])))
    .on("mousemove", moveTooltip)
    .on("mouseleave", hideTooltip);

  cell.append("text")
    .attr("x", cellWidth / 2)
    .attr("y", cellHeight / 2)
    .attr("dy", "0.35em")
    .attr("text-anchor", "middle")
    .attr("class", "value-label")
    .attr("fill", (d) => Number(d.row[d.key] || 0) > 100000 ? "#fff" : "#51615b")
    .text((d) => {
      const value = Number(d.row[d.key] || 0);
      return value > 0 ? formatCompact(value) : "";
    });
}

function renderAnnotationTerms(data) {
  const rows = (data.annotation_terms[state.annotationDomain] || [])
    .slice(0, 24)
    .map((item) => ({
      key: item.term_id,
      label: item.label === item.term_id ? item.term_id : `${item.label}`,
      value: Number(item.global_count || 0),
      color: categorical(item.ontology_prefix),
      tooltip: tooltipHtml(item.label, [
        ["Term", item.term_id],
        ["Ontology", item.ontology_prefix],
        ["Count", formatNumber(item.global_count)]
      ])
    }));

  drawHorizontalBars("#annotation-terms", rows);
}

function renderTableInventory(data) {
  const rows = data.table_inventory
    .slice(0, 14)
    .map((table) => ({
      key: table.relation,
      label: table.relation,
      value: Number(table.total_size_bytes || 0),
      color: table.kind === "materialized view" ? palette.green : palette.blue,
      tooltip: tooltipHtml(table.relation, [
        ["Kind", table.kind],
        ["Rows", formatNumber(table.estimated_rows)],
        ["Size", formatBytes(table.total_size_bytes)]
      ])
    }));

  drawHorizontalBars("#table-inventory", rows, {
    valueFormatter: (value) => formatBytes(value)
  });
}

function renderAll() {
  const data = state.data;
  if (!data) return;

  renderKpis(data);
  renderResourceLandscape(data);
  renderRelationFlow(data);
  renderResourceBars(data);
  renderEntityFacet(data);
  renderRelationFacet(data);
  renderResourceMatrix(data);
  renderAnnotationTerms(data);
  renderTableInventory(data);
}

function setActiveButton(selector, activeValue, attribute) {
  d3.selectAll(selector).classed("active", function () {
    return this.dataset[attribute] === activeValue;
  });
}

function bindControls() {
  d3.select("#resource-metric").on("change", (event) => {
    state.resourceMetric = event.target.value;
    renderResourceBars(state.data);
  });

  d3.selectAll("[data-entity-facet]").on("click", (event) => {
    state.entityFacet = event.currentTarget.dataset.entityFacet;
    setActiveButton("[data-entity-facet]", state.entityFacet, "entityFacet");
    renderEntityFacet(state.data);
  });

  d3.selectAll("[data-relation-facet]").on("click", (event) => {
    state.relationFacet = event.currentTarget.dataset.relationFacet;
    setActiveButton("[data-relation-facet]", state.relationFacet, "relationFacet");
    renderRelationFacet(state.data);
  });

  d3.selectAll("[data-annotation-domain]").on("click", (event) => {
    state.annotationDomain = event.currentTarget.dataset.annotationDomain;
    setActiveButton("[data-annotation-domain]", state.annotationDomain, "annotationDomain");
    renderAnnotationTerms(state.data);
  });

  let resizeFrame = null;
  window.addEventListener("resize", () => {
    cancelAnimationFrame(resizeFrame);
    resizeFrame = requestAnimationFrame(renderAll);
  });
}

async function init() {
  if (!d3) {
    throw new Error("D3 did not load.");
  }

  bindControls();
  state.data = await d3.json("./data/figures.json");
  renderAll();
}

init().catch((error) => {
  console.error(error);
  const shell = document.querySelector(".shell");
  const message = document.createElement("div");
  message.className = "error-state";
  message.textContent = "Unable to load the generated figure data.";
  shell?.append(message);
});
