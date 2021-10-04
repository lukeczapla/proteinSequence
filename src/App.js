import logo from './logo.svg';
import {useState, useEffect, useRef} from 'react';
import './App.css';
import Card from 'react-bootstrap/Card';
import Row from 'react-bootstrap/Row';
import Col from 'react-bootstrap/Col';
import Form from 'react-bootstrap/Form';
import Container from 'react-bootstrap/Container';
import Button from 'react-bootstrap/Button';
import Table from 'react-bootstrap/Table';
import 'bootstrap/dist/css/bootstrap.min.css';
import * as d3 from 'd3';

const App = (props) => {

    const svgRef = useRef(null);

    const [state, setState] = useState({
        epitopeCentral: [[], [], [], []],
        epitopeCentralNorm: [[], [], [], []],
        epitopeSliding: [[], [], [], []],
        epitopeSlidingNorm: [[], [], [], []],
        antigenSequenceCentral: [[], [], [], []],
        antigenSequenceCentralNorm: [[], [], [], []],
        antigenSequenceSliding: [[], [], [], []],
        antigenSequenceSlidingNorm: [[], [], [], []],
        total: 0,
        dataSet: null,
        dataRow: null,
        colors: null,
        modelSelected: "none",
        plotSelected: "aaTypes",
        plot: null
    });

    const dataNames = [[["A","C","D","E","F","G","H","I","K","L","M","N","P","Q","R","S","T","V","W","Y"], "Alanine, Cytosine, Aspartic Acid, Glutamic Acid, Phenylalanine, Glycine, Histidine, Isoleucine, " +
    "Lysine, Leucine, Methionine, Asparagine, Proline, Glutamine, Arginine, Serine, Threonine, Valine, Tryptophan, Tyrosine"], [["A","N","P","C"], "Anionic (-), Non-Polar, Polar (Uncharged), Cationic (+)"]]


    useEffect(() => {
        let value = "composition.txt";
        let readEpitope = [[], [], [], []];
        fetch(value).then(data => data.text()).then(result => {
            let lines = result.split("\n");
            lines.forEach((line, i) => {
                //console.log(i + ": " + line);
                if (i < 4) {
                    let vals = line.split(" ");
                    for (let j = 0; j < 21; j++) {
                        readEpitope[i].push(parseFloat(vals[j]));
                    }
                    if (i === 3) {
                        state.epitopeCentral = JSON.parse(JSON.stringify(readEpitope));
                        readEpitope = [[], [], [], []];
                    }
                }
                if (i >= 5 && i <= 8) {
                    let vals = line.split(" ");
                    for (let j = 0; j < 21; j++) {
                        readEpitope[i-5].push(parseFloat(vals[j]));
                    }
                    if (i === 8) {
                        state.epitopeCentralNorm = JSON.parse(JSON.stringify(readEpitope));
                        readEpitope = [[], [], [], []];
                    }
                }
                if (i >= 10 && i <= 13) {
                    let vals = line.split(" ");
                    for (let j = 0; j < 21; j++) {
                        readEpitope[i-10].push(parseFloat(vals[j]));
                    }
                    if (i === 13) {
                        state.epitopeSliding = JSON.parse(JSON.stringify(readEpitope));
                        readEpitope = [[], [], [], []];
                    }
                }
                if (i >= 15 && i <= 18) {
                    let vals = line.split(" ");
                    for (let j = 0; j < 21; j++) {
                        readEpitope[i-15].push(parseFloat(vals[j]));
                    }
                    if (i === 18) {
                        state.epitopeSlidingNorm = JSON.parse(JSON.stringify(readEpitope));
                        readEpitope = [[], [], [], []];
                    }
                }

                if (i >= 22 && i <= 25) {
                    let vals = line.split(" ");
                    for (let j = 0; j < 21; j++) {
                        readEpitope[i-22].push(parseFloat(vals[j]));
                    }
                    if (i === 25) {
                        state.antigenSequenceCentral = JSON.parse(JSON.stringify(readEpitope));
                        readEpitope = [[], [], [], []];
                    }
                }
                if (i >= 27 && i <= 30) {
                    let vals = line.split(" ");
                    for (let j = 0; j < 21; j++) {
                        readEpitope[i-27].push(parseFloat(vals[j]));
                    }
                    if (i === 30) {
                        state.antigenSequenceCentralNorm = JSON.parse(JSON.stringify(readEpitope));
                        readEpitope = [[], [], [], []];
                    }
                }
                if (i >= 32 && i <= 35) {
                    let vals = line.split(" ");
                    for (let j = 0; j < 21; j++) {
                        readEpitope[i-32].push(parseFloat(vals[j]));
                    }
                    if (i === 35) {
                        state.antigenSequenceSliding = JSON.parse(JSON.stringify(readEpitope));
                        readEpitope = [[], [], [], []];
                    }
                }
                if (i >= 37 && i <= 40) {
                    let vals = line.split(" ");
                    for (let j = 0; j < 21; j++) {
                        readEpitope[i-37].push(parseFloat(vals[j]));
                    }
                    if (i === 40) {
                        state.antigenSequenceSlidingNorm = JSON.parse(JSON.stringify(readEpitope));
                        readEpitope = [[], [], [], []];
                    }
                }
            });
        });
    }, []);


    useEffect(() => {
        if (state.modelSelect !== "none") process(state.modelSelected);
    }, [state.modelSelected, state.plotSelected]);

    const aminoAcids = "ACDEFGHIKLMNPQRSTVWY";

    const aminoAcidMap = {anionic: "DE", nonpolar: "AVLGIMWFP", polar: "CNQSTY", cationic: "KRH"};
    const groupList = ["anionic", "non-polar", "polar", "cationic"];
    const codeMap = {0: "anionic", 1: "non-polar", 2: "polar", 3: "cationic"};

    const grouper = (letter) => {
        if (aminoAcidMap.nonpolar.includes(letter)) return [1, "Non-Polar"];
        if (aminoAcidMap.polar.includes(letter)) return [2, "Polar"];
        if (aminoAcidMap.anionic.includes(letter)) return [0, "Anionic"];
        if (aminoAcidMap.cationic.includes(letter)) return [3, "Cationic"];
    }



    const prepData = (m, ds, group = false) => {
        //console.log(m);
        //console.log(group);
        let colors = d3.schemeSpectral[group ? 4: 10];
        let dataSet = ds;
        let total2 = 0;
        let dataRow = null;
        if (group) dataRow = 1;
        else dataRow = 0;
        let data = [];
        let ds2 = [[], [], [], []];
        for (let v = 0; v < 4; v++) {
            let obj = {};
            if (v === 0) obj["name"] = "1st Amino Acid";
            if (v === 1) obj["name"] = "2nd Amino Acid";
            if (v === 2) obj["name"] = "3rd Amino Acid";
            if (v === 3) obj["name"] = "4th Amino Acid";
            total2 = 0;
            if (!group) {
                for (let w = 0; w < 20; w++) {
                    obj[aminoAcids[w]] = m[v][w];
                }
                obj["total"] = m[v][20];
            } else {
                let data = [0, 0, 0, 0];

                for (let w = 0; w < 20; w++) {
                    data[grouper(aminoAcids[w])[0]] += m[v][w];
                }
                total2 = 0;
                for (let w = 0; w < 4; w++) {
                    total2 += data[w];
                    obj[codeMap[w]] = data[w];
                }
                for (let w = 0; w < 4; w++) {
                    ds2[v].push(data[w]/total2);
                }
                obj["total"] = total2;
            }
            data.push(obj);
        }
        let arr = ["name"];
        if (group) setState({...state, colors: colors, dataSet: ds2, dataRow: dataRow, total: total2});
        else setState({...state, colors: colors, dataSet: dataSet, dataRow: dataRow, total: m[0][20]});

        if (!group) for (let v = 0; v < 20; v++) {
            arr.push(aminoAcids[v]);
        } else for (let v = 0; v < 4; v++) {
            arr.push(groupList[v]);
        }
        data["columns"] = arr;
        let margin = {top: 30, right: 10, bottom: 0, left: 30};
        let width = 600;
        let height = data.length * 25 + margin.top + margin.bottom
        let series = d3.stack()
            .keys(data.columns.slice(1))
            .offset(d3.stackOffsetExpand)
            (data)
            .map(d => (d.forEach(v => v.key = d.key), d));
        let x = d3.scaleLinear()
            .range([margin.left, width - margin.right]);
        let y = d3.scaleBand()
            .domain(data.map(d => d.name))
            .range([margin.top, height - margin.bottom])
            .padding(0.08);
        let xAxis = g => g
            .attr("transform", `translate(0,${margin.top})`)
            .call(d3.axisTop(x).ticks(width / 100, "%"))
            .call(g => g.selectAll(".domain").remove())
        let yAxis = g => g
            .attr("transform", `translate(${margin.left},0)`)
            .call(d3.axisLeft(y).tickSizeOuter(0))
            .call(g => g.selectAll(".domain").remove());
        //console.log("DATA");
        //console.log(data);
        //console.log("SERIES");
        //console.log(series);
        let formatValue = x => isNaN(x) ? "N/A" : x.toLocaleString("en");
        let formatPercent = d3.format(".1%");
        let color = d3.scaleOrdinal()
            .domain(series.map(d => d.key))
            .range(d3.schemeSpectral[group ? 4: 10])
            .unknown("#ccc");
        //let key = legend({title: "Age (years)", color, tickSize: 0});
        let svgE = d3.select(svgRef.current);
        svgE.selectAll("*").remove();
        let svg = svgE
            .attr("viewBox", [0, 0, width, height])
            .style("overflow", "visible");
        svg.append("g")
            .selectAll("g")
            .data(series)
            .enter().append("g")
            .attr("fill", d => color(d.key))
            .selectAll("rect")
            .data(d => d)
            .join("rect")
            .attr("x", d => x(d[0]))
            .attr("y", (d, i) => y(d.data.name))
            .attr("width", d => x(d[1]) - x(d[0]))
            .attr("height", y.bandwidth())
            .append("title")
                    .text(d => `${d.data.name} ${d.key}
                ${formatPercent(d[1] - d[0])} (${formatValue(d.data[d.key])})`);

        svg.append("g")
            .call(xAxis);

        svg.append("g")
            .call(yAxis);
        console.log(d3.schemeSpectral[group ? 4: 10]);
        return svg.node();

    }

    const changeEvent = (e) => {
        if (e.target.name === "modelSelected") {
            setState({...state, modelSelected: e.target.value});
        }
        if (e.target.name === "plotSelected") {
            setState({...state, plotSelected: e.target.value});
        }
    }

    const process = (modelSelected) => {
        switch (modelSelected) {
            case "none": return;
            case "epitopeCentral": return prepData(state.epitopeCentral, state.epitopeCentralNorm, state.plotSelected === "aaCategories");
            case "epitopeSliding": return prepData(state.epitopeSliding, state.epitopeSlidingNorm, state.plotSelected === "aaCategories");
            case "antigenSequenceCentral": return prepData(state.antigenSequenceCentral, state.antigenSequenceCentralNorm, state.plotSelected === "aaCategories");
            case "antigenSequenceSliding": return prepData(state.antigenSequenceSliding, state.antigenSequenceSlidingNorm, state.plotSelected === "aaCategories");
            default: return;
        }
    }


  return (
    <div className="App">
      <header className="App-header">
          <Row>
            <Card style={{fontFamily: "Segoe UI", color:"lightblue", backgroundColor:"#282c34"}}>Welcome to The Language of Biology 1: Protein Sequence Magic</Card>
          </Row>
      </header>
      <Container>
        <Row className="justify-content-md-center">
          <Col>
              <Form.Select id="mselect" name="modelSelected" onChange={(e) => changeEvent(e)} value={state.modelSelected}>
                  <option value="none" key={0}>--Select a model below--</option>
                  <option value="epitopeCentral" key={1}>Epitope Central amino acids</option>
                  <option value="epitopeSliding" key={2}>Epitope Sliding Window amino acids</option>
                  <option value="antigenSequenceCentral" key={3}>Antigen Sequence Central amino acids</option>
                  <option value="antigenSequenceSliding" key={4}>Antigen Sequence Sliding Window amino acids</option>
              </Form.Select>
          </Col>
          <Col>
              <Form.Select id="pselect" name="plotSelected" onChange={(e) => changeEvent(e)} value={state.plotSelected}>
                  <option value="aaTypes" key={1}>Amino Acids</option>
                  <option value="aaCategories" key={2}>Amino Acid Categories</option>
              </Form.Select>
          </Col>
        </Row>
        <Row className="justify-content-md-center">
            <Button onClick={() => process(state.modelSelected)}>Process Results</Button>
        </Row>
      </Container>
        <svg ref={svgRef} width={900} height={400} />
      <Container>
          {svgRef && state.dataRow !== null && state.dataSet !== null &&
          <Table striped bordered hover>
          <thead>
          <tr>
              <th>{state.dataRow === 0 ? "Letter" : "Code"}</th>
              <th>Name</th>
              <th>Percentages in AA 1</th>
              <th>Percentages in AA 2</th>
              <th>Percentages in AA 3</th>
              <th>Percentages in AA 4</th>
              <th>Total Samples</th>
          </tr>
          </thead>
          <tbody>
          {state.dataRow !== null && dataNames[state.dataRow][0].map((letter, i) => (
              <tr>
                  <td style={{backgroundColor: state.colors[i % 10], color: "white"}}><b>{letter}</b></td>
                  <td><b>{dataNames[state.dataRow][1].split(", ")[i]}</b></td>
                  <td>{parseFloat(100.0*state.dataSet[0][i]).toFixed(1) + "%"}</td>
                  <td>{parseFloat(100.0*state.dataSet[1][i]).toFixed(1) + "%"}</td>
                  <td>{parseFloat(100.0*state.dataSet[2][i]).toFixed(1) + "%"}</td>
                  <td>{parseFloat(100.0*state.dataSet[3][i]).toFixed(1) + "%"}</td>
                  {i === 0 && <td>{i === 0 ? <b>{state.total}</b>: "..."}</td>}
              </tr>
          ))}
          </tbody>
          </Table>}
      </Container>
    </div>
  );
}

export default App;
