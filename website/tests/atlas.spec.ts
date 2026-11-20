import { expect, test, type Page } from "@playwright/test";

const origin = "http://127.0.0.1:4173";
const results = (page: Page) => page.locator(".result-card");
const select = (page: Page, id: string) => page.locator(".region-index").getByRole("button", { name: new RegExp(" " + id + "，") });

test.beforeEach(async ({page}) => {
  await page.route("**/*", route => route.request().url().startsWith(origin) ? route.continue() : route.abort());
  await page.goto("/");
  await expect(results(page).getByRole("heading", { name: "丘脑前背核" })).toBeVisible();
  await expect(results(page).locator(".big-effect")).toContainText("+0.908");
});

test("mouse atlas switches decoded regions through SVG, keyboard and synchronized views", async ({page}) => {
  const map = page.locator(".brain-image");
  const target = map.getByRole("button", { name: "海马 CA1 区 CA1", exact: true });
  await expect(target.locator("title")).toHaveText(/Hippocampal field CA1/);
  await target.focus();
  await expect(page.locator(".map-focus")).toContainText("海马 CA1 区");
  await page.keyboard.press("Enter");
  await expect(results(page).locator(".big-effect")).toContainText("+0.837");
  await select(page, "PoSub").click();
  await expect(results(page).getByRole("heading", { name: "后下托" })).toBeVisible();
  await expect(results(page).locator(".big-effect")).toContainText("+1.319");
  await page.getByRole("button", { name: "冠状面", exact: true }).click();
  await expect(map.getByRole("button", { name: "后下托 PoSub", exact: true })).toHaveAttribute("aria-pressed", "true");
  const pag = map.getByRole("button", { name: "导水管周围灰质 PAG", exact: true });
  const [cx, cy] = (await pag.getAttribute("data-center"))!.split(",").map(Number);
  const point = await pag.evaluate((el, [x, y]) => {
    const matrix = (el as unknown as SVGPathElement).getScreenCTM()!;
    const p = new DOMPoint(x, y).matrixTransform(matrix);
    return {x: p.x, y: p.y};
  }, [cx, cy]);
  await page.mouse.click(point.x, point.y);
  await expect(results(page).getByRole("heading", { name: "导水管周围灰质" })).toBeVisible();
  await expect(results(page).locator(".big-effect")).toContainText("+0.257");
});

test("PoSub source selection and contrasts preserve denominator and run lineage", async ({page}) => {
  await select(page, "PoSub").click();
  const panel = results(page);
  await expect(panel.getByLabel("数据来源")).toHaveValue("000939");
  await expect(panel.locator(".big-effect")).toContainText("+1.319");
  await expect(panel.locator(".chart-intro")).toContainText("29 动物");
  await panel.getByRole("button", { name: "NREM / REM", exact: true }).click();
  await expect(panel.locator(".big-effect")).toContainText("−1.174");
  await expect(panel.locator(".chart-intro")).toContainText("25 动物");
  await expect(page.getByRole("link", { name: "原始论文", exact: true })).toHaveAttribute("href", "https://doi.org/10.1038/s41593-024-01588-5");
  await expect(page.getByRole("link", { name: "原始数据", exact: true })).toHaveAttribute("href", /dandi\.000939\/0\.261016\.2310/);
  await panel.getByLabel("数据来源").selectOption("000056");
  await expect(panel.locator(".big-effect")).toContainText("−0.363");
  await expect(panel.locator(".chart-intro")).toContainText("3 动物");
  await panel.getByLabel("选择动物").selectOption("Mouse28");
  await expect(panel.locator(".session-drill")).toContainText("Mouse28");
  await panel.locator(".session-drill summary").first().click();
  await expect(panel.locator(".fit-list").first()).toContainText("bits/s");
});

test("missing REM remains missing and mouse and rat data are separate", async ({page}) => {
  await select(page, "V1").click();
  const panel = results(page);
  await expect(panel.locator(".chart-intro")).toContainText("13 会话");
  await expect(panel.locator(".region-metrics")).toContainText("未核实");
  await panel.getByRole("button", { name: "Wake / REM" }).click();
  await expect(panel).toContainText("无已核验的状态对比");
  await panel.getByRole("button", { name: "三态", exact: true }).click();
  await expect(panel).toContainText("缺少可用状态：REM");
  await expect(panel.locator(".chart-replicate circle")).toHaveCount(26);
  await page.getByRole("button", { name: "大鼠", exact: true }).click();
  await expect(page.locator(".brain-image")).toHaveCount(0);
  await expect(page.getByRole("heading", { name: "大鼠数据探索" })).toBeVisible();
  await expect(panel.getByLabel("数据来源")).toHaveValue("000041");
  await expect(page.getByRole("link", { name: "原始论文", exact: true })).toHaveAttribute("href", "https://doi.org/10.1016/j.neuron.2016.03.036");
  await expect(select(page, "V1")).toHaveCount(0);
  await page.getByRole("button", { name: "小鼠", exact: true }).click();
  await expect(page.locator(".brain-image")).toBeVisible();
  await select(page, "RSC").click();
  await expect(panel).toContainText("现有记录未获得可展示的状态估计");
});

test("model controls and model adequacy use the saved paired diagnostics", async ({page}) => {
  const qc = page.locator("#quality");
  await expect(qc.getByLabel("分析规格")).toHaveValue("AR10");
  await expect(qc.locator(".method-comparison")).toContainText("+0.908");
  await expect(qc.locator(".method-comparison")).toContainText("+0.640");
  await expect(qc.locator(".qc-notes")).toContainText("模型充分性提示");
  await expect(qc.locator(".qc-notes")).toContainText("残差");
  await qc.getByLabel("分析规格").selectOption("rate_thinning");
  await expect(qc.locator(".method-comparison")).toContainText("配平时间窗的参考");
  await expect(qc).toContainText("对照采用各动物的重复中位差");
  await expect(results(page).locator(".big-effect")).toContainText("+0.908");
});

test("responsive page has no horizontal clipping at 390 pixels and serves all data locally", async ({page}) => {
  const requests: string[] = [];
  page.on("request", req => requests.push(req.url()));
  await page.setViewportSize({width: 390, height: 844});
  await page.reload();
  await expect(page.getByRole("heading", {name: "睡眠临界性脑图谱"})).toBeVisible();
  await select(page, "PAG").focus();
  await page.keyboard.press("Space");
  await expect(results(page).locator(".big-effect")).toContainText("+0.257");
  await expect(page.locator(".brain-image")).toBeVisible();
  expect(await page.evaluate(() => document.documentElement.scrollWidth)).toBe(390);
  expect(requests.length).toBeGreaterThan(1);
  expect(requests.every(x => x.startsWith(origin) || x.startsWith("data:"))).toBeTruthy();
});

test("audit sensitivity is visible and the single-state rat source keeps AR methods distinct", async ({page}) => {
  await select(page, "V1").click();
  await expect(results(page).locator(".big-effect")).toContainText("+1.601");
  const qc = page.locator("#quality");
  await qc.getByText(/针对性审计/).click();
  await expect(qc).toContainText("AR(10) +0.939 bits/s");
  await page.getByRole("button", { name: "大鼠", exact: true }).click();
  await select(page, "mPFC").click();
  await results(page).getByLabel("数据来源").selectOption("hc-24");
  await expect(results(page)).toContainText("无已核验的状态对比");
  await expect(qc).toContainText("NREM d₂ · 同记录对照");
  await expect(qc.getByLabel("分析规格")).toHaveValue("AR10");
  await expect(qc.locator(".method-comparison").getByText("AR(10)", {exact:true})).toBeVisible();
  await page.getByRole("button", {name:"三态",exact:true}).click();
  await expect(results(page).locator(".chart-replicate circle")).toHaveCount(1);
});
