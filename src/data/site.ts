export const site = {
  name: 'Isaque Dutra',
  role: 'PhD student',
  department: 'Department of Physics',
  institution: 'Yale University',
  field: 'Astrophysics & Cosmology',
  description:
    'Research interests and publications in astrophysics and cosmology: dark matter, gravitational lensing, and the expansion of the Universe.',
  lead: 'I study dark matter and the expansion of the Universe through gravitational lensing and AGN variability.',
  email: 'isaque.dutra@yale.edu',
  url: 'https://dutra.cc',
  links: [
    {
      label: 'arXiv',
      href: 'https://arxiv.org/search/?query=Dutra%2C+Isaque&searchtype=author&abstracts=show&order=-announced_date_first&size=50',
    },
    {
      label: 'ADS',
      href: 'https://ui.adsabs.harvard.edu/search/q=author%3A%22Dutra%2C%20Isaque%22&sort=date%20desc%2C%20bibcode%20desc',
    },
    { label: 'ORCID', href: 'https://orcid.org/0000-0001-7040-4930' },
  ],
  nav: [
    { label: 'About', href: '/#about' },
    { label: 'Research', href: '/#research' },
    { label: 'Publications', href: '/#publications' },
    { label: 'Projects', href: '/projects/' },
    { label: 'Contact', href: '/#contact' },
  ],
};

export const about = [
  'I use galaxy clusters, the most massive gravitationally bound structures in the universe, as laboratories for studying dark matter. By comparing the dark matter distributions inferred from gravitational lensing with cosmological simulations, I investigate whether discrepancies between observations and predictions reflect limitations in our models or point toward alternate dark matter physics.',
  'With collaborators, I have also used the variability of active galactic nuclei (AGN) as a cosmological distance probe to investigate the cosmic expansion and the possible evolution of dark energy.',
];

export const research = [
  {
    topic: 'Dark matter',
    question: 'Why do cluster galaxies lens more strongly than cold dark matter predicts?',
    answer:
      'Observed galaxy clusters produce nearly ten times more galaxy–galaxy strong lensing than ΛCDM simulations predict, and I use this discrepancy to probe the nature of dark matter. By measuring how strongly each member galaxy lenses, both on its own and embedded in the cluster, we have shown that the gap has two sources: simulated galaxies are not centrally concentrated enough, and too few of them sit where the cluster’s lensing field is strongest. We investigated whether self-interacting dark matter (SIDM), whose subhalos can undergo gravothermal core collapse and grow steep, dense centres, could supply the missing concentration. A contraction that mimics this core collapse recovers the first effect but not the second, pointing to differences in both the inner structure and the locations of cluster galaxies.',
  },
  {
    topic: 'Distance probes',
    question: 'Can AGN variability measure cosmic distances?',
    answer:
      'AGN variability may offer a way to measure distances across cosmic history. With collaborators, I investigate how the relationship between variability and intrinsic luminosity can be used for cosmology. The goal is to develop a reliable distance probe that reaches high redshift and to understand the assumptions and uncertainties that shape its conclusions.',
  },
  {
    topic: 'Dark energy',
    question: 'Does dark energy evolve over time?',
    answer:
      'I use the distance ratios encoded in gravitational lensing and the distances inferred from AGN variability to investigate how the Universe has expanded over time. These approaches provide complementary ways to test whether the expansion history is consistent with a cosmological constant or suggests that dark energy has evolved since the Big Bang.',
  },
];
